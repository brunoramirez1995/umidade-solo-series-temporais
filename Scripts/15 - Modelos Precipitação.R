# 15 - Modelos AR com precipitação e interação

# Pacotes e configurações

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)


source("Scripts/apoio_modelagem.R", encoding = "UTF-8")

arquivo_dados <- "Dados/nasa_power_diario.csv"
pasta_resultados <- "resultados/modelagem/precipitacao"
pasta_figuras <- "figuras/principais"
chaves_serie <- c("localidade_id", "localidade", "clima", "variavel")
chaves_modelo <- c("localidade_id", "variavel", "modelo")


inicio_grafico <- as.Date("2024-10-01")
fim_grafico <- as.Date("2024-12-31")
localidades_grafico <- c("Belém", "Brasília", "Petrolina")
cores <- c("AR" = "#009E73", "AR + P" = "#2166AC",
           "AR + P + interação" = "#D55E00", "Observado" = "grey25")


# Configurações e previsões já obtidas nos scripts 12 e 13

referencia_ar <- read_csv(
  "resultados/modelagem/autorregressivos/selecionados_validacao.csv",
  show_col_types = FALSE
) |>
  select(all_of(chaves_serie), modelo_AR = modelo, p, K, tendencia)

previsoes_ar <- read_csv(
  "resultados/modelagem/autorregressivos/previsoes_validacao.csv",
  show_col_types = FALSE
) |>
  semi_join(referencia_ar,
            by = c("localidade_id", "variavel", "modelo" = "modelo_AR")) |>
  mutate(referencia = "AR")

persistencia <- read_csv(
  "resultados/modelagem/referencias/previsoes_validacao.csv",
  show_col_types = FALSE
) |>
  filter(modelo == "persistencia") |>
  mutate(referencia = "persistencia")

referencias <- bind_rows(previsoes_ar, persistencia)


# Leitura e construção das entradas

dados <- read_csv(
  arquivo_dados,
  col_types = cols_only(
    localidade_id = col_character(), localidade = col_character(),
    clima = col_character(), data = col_date(),
    GWETTOP = col_double(), GWETROOT = col_double(),
    PRECTOTCORR = col_double()
  ),
  show_col_types = FALSE
) |>
  filter(data >= inicio_treinamento, data <= fim_validacao)

# A precipitação é a nova entrada

dados <- dados |>
  pivot_longer(all_of(alvos), names_to = "variavel", values_to = "observado") |>
  arrange(localidade_id, variavel, data) |>
  group_by(localidade_id, variavel) |>
  mutate(
    lag_1 = lag(observado, 1),
    lag_2 = lag(observado, 2),
    lag_3 = lag(observado, 3),
    precipitacao_lag_1 = lag(PRECTOTCORR, 1),
    interacao = lag_1 * precipitacao_lag_1
  ) |>
  ungroup() |>
  mutate(tempo = as.numeric(data - inicio_treinamento) / 365.25)

for (k in 1:3) {
  dados[[paste0("sen_", k)]] <- sin(2 * pi * k * dados$tempo)
  dados[[paste0("cos_", k)]] <- cos(2 * pi * k * dados$tempo)
}

# Cada linha corresponde ao dia previsto, t+1.
# lag_1 e precipitacao_lag_1 são Y_t e P_t; a interação é Y_t * P_t.
# Os harmônicos usam a data prevista. A chuva desse dia não entra no modelo.
# A precipitação permanece em mm/dia, sem transformação ou classificação em faixas.


# Duas extensões do AR selecionado em cada série

configuracoes <- referencia_ar |>
  crossing(com_interacao = c(FALSE, TRUE)) |>
  mutate(
    familia = if_else(com_interacao, "AR + P + interação", "AR + P"),
    modelo = paste0(modelo_AR, if_else(com_interacao, "_P_INTERACAO", "_P")),
    n_parametros = 3L + p + 2L * K + as.integer(tendencia) + as.integer(com_interacao)
  )

# p, K e a presença de tendência vêm do script 13. Não há uma nova grade.
# Todos os coeficientes serão reestimados, inclusive os lags e os harmônicos.
# São 72 ajustes: duas extensões, 12 séries e três janelas anuais.
# A contagem de parâmetros inclui a variância dos erros.


# Ajustes anuais e previsões de um dia à frente

lista_previsoes <- list()
lista_ajustes <- list()
lista_coeficientes <- list()
falhas <- tibble(
  localidade_id = character(), variavel = character(), modelo = character(),
  ano = integer(), motivo = character()
)

for (i in seq_len(nrow(configuracoes))) {
  config <- configuracoes[i, ]
  message(config$localidade, " - ", config$variavel, " - ", config$familia)

  serie <- dados |>
    filter(localidade_id == config$localidade_id, variavel == config$variavel)
  entradas <- c(
    paste0("lag_", seq_len(config$p)),
    paste0("sen_", seq_len(config$K)),
    paste0("cos_", seq_len(config$K)),
    if (config$tendencia) "tempo",
    "precipitacao_lag_1",
    if (config$com_interacao) "interacao"
  )
  formula_modelo <- reformulate(entradas, response = "observado")

  for (j in seq_len(nrow(calendario_validacao))) {
    janela <- calendario_validacao[j, ]
    treino <- serie |>
      filter(data >= janela$inicio_ajuste, data <= janela$fim_ajuste) |>
      drop_na(observado, all_of(entradas))
    avaliacao <- serie |>
      filter(data >= janela$inicio_avaliacao, data <= janela$fim_avaliacao)

    resultado <- tryCatch({
      ajuste <- lm(formula_modelo, data = treino)
      if (any(!is.finite(coef(ajuste)))) stop("Coeficientes não identificados.")

      previsto <- as.numeric(predict(ajuste, newdata = avaliacao))
      if (any(!is.finite(previsto))) stop("Previsões não finitas.")

      list(ajuste = ajuste, previsoes = avaliacao |>
             transmute(localidade_id, localidade, clima, variavel,
                       data_origem = data - 1, data_previsao = data,
                       ano = janela$ano, familia = config$familia, modelo = config$modelo,
                       observado, previsto = .env$previsto, erro = observado - previsto))
    }, error = function(e) e)

    if (inherits(resultado, "error")) {
      falhas <- bind_rows(falhas, tibble(
        localidade_id = config$localidade_id, variavel = config$variavel,
        modelo = config$modelo, ano = janela$ano, motivo = conditionMessage(resultado)
      ))
    } else {
      lista_previsoes[[length(lista_previsoes) + 1L]] <- resultado$previsoes
      lista_ajustes[[length(lista_ajustes) + 1L]] <- tibble(
        localidade_id = config$localidade_id, variavel = config$variavel,
        familia = config$familia, modelo = config$modelo,
        p = config$p, K = config$K, tendencia = config$tendencia,
        com_interacao = config$com_interacao, ano = janela$ano,
        inicio_ajuste = janela$inicio_ajuste, fim_ajuste = janela$fim_ajuste,
        dados_treino = list(select(treino, data, observado)),
        ajuste = list(resultado$ajuste)
      )

      coeficientes <- coef(resultado$ajuste)
      lista_coeficientes[[length(lista_coeficientes) + 1L]] <- tibble(
        localidade_id = config$localidade_id, localidade = config$localidade,
        variavel = config$variavel, familia = config$familia, modelo = config$modelo,
        ano = janela$ano, beta = unname(coeficientes["precipitacao_lag_1"]),
        gamma = if (config$com_interacao) unname(coeficientes["interacao"]) else 0
      )
    }
  }
}

previsoes <- bind_rows(lista_previsoes)
coeficientes_precipitacao <- bind_rows(lista_coeficientes)

# Os coeficientes ficam fixos no ano. Os valores observados atualizam as entradas.

# Os ajustes são guardados aqui para os diagnósticos


# Métricas gerais e anuais, comparadas ao AR e à persistência

resumir_previsoes <- function(grupos) {
  metricas <- previsoes |>
    group_by(across(all_of(c(grupos, "familia", "modelo")))) |>
    summarise(calcular_metricas(observado, previsto),
              n_fora = sum(previsto < 0 | previsto > 1),
              percentual_fora = 100 * n_fora / n,
              excesso_maximo = max(pmax(-previsto, previsto - 1, 0)), .groups = "drop")

  comparadores <- referencias |>
    group_by(across(all_of(c(grupos, "referencia")))) |>
    summarise(RMSE = sqrt(mean(erro^2)), MAE = mean(abs(erro)), .groups = "drop") |>
    pivot_wider(names_from = referencia, values_from = c(RMSE, MAE))

  metricas |>
    left_join(comparadores, by = grupos) |>
    mutate(ganho_RMSE_AR = 100 * (1 - RMSE / RMSE_AR),
           ganho_MAE_AR = 100 * (1 - MAE / MAE_AR),
           ganho_RMSE_persistencia = 100 * (1 - RMSE / RMSE_persistencia),
           ganho_MAE_persistencia = 100 * (1 - MAE / MAE_persistencia))
}


metricas <- resumir_previsoes(chaves_serie) |>
  left_join(select(configuracoes, all_of(chaves_modelo), p, K, tendencia,
                   com_interacao, n_parametros), by = chaves_modelo) |>
  mutate(completo = n == as.integer(fim_validacao - inicio_validacao) + 1L,
         across(starts_with("ganho_"), ~ if_else(completo, .x, NA_real_)))


metricas_anuais <- resumir_previsoes(c(chaves_serie, "ano"))

# Cada família tem uma configuração por série. Somente as completas são mantidas.


selecionados <- metricas |> filter(completo)
ajustes_selecionados <- bind_rows(lista_ajustes) |>
  semi_join(selecionados, by = chaves_modelo)

print(metricas |>
        select(localidade, variavel, familia, RMSE, MAE, ganho_RMSE_AR,
               ganho_MAE_AR, ganho_RMSE_persistencia, vies, n_fora, completo),
      n = Inf, width = Inf)


# A precipitação reduziu o RMSE de GWETROOT nas seis localidades.
# Em GWETTOP, acrescentar apenas a chuva pouco alterou o desempenho.
# A interação trouxe o maior ganho em Petrolina, seguida por Montes Claros e Brasília.
# Frente ao modelo com chuva sem interação, ela melhorou o RMSE em 11 séries
# e o MAE em oito. Portanto, seu ganho não foi uniforme entre as métricas.



print(metricas_anuais |>
        select(localidade, variavel, familia, ano, RMSE, MAE,
               ganho_RMSE_AR, ganho_MAE_AR, vies, n_fora), n = Inf, width = Inf)

# O ganho da chuva em GWETROOT apareceu nos três anos e nas seis localidades.
# Em Petrolina, a interação melhorou o RMSE de GWETTOP nos três anos.
# Houve exceções anuais na superfície: Manaus em 2022 e Uberlândia em 2023.


print(coeficientes_precipitacao, n = Inf, width = Inf)

print(falhas, n = Inf, width = Inf)

# Ganho positivo indica redução do erro. RMSE e MAE podem levar a conclusões distintas.
# Para avaliar a interação, compare AR + P + interação com AR + P na mesma série.
# Confira também os três anos: uma melhora agregada pode depender de um único ano.
# O viés é observado - previsto. Valores negativos indicam superestimação média.

# No modelo com interação, a contribuição por mm de chuva é beta + gamma * Y_t.
# gamma negativo indica que essa contribuição diminui conforme a umidade aumenta;
# não significa, sozinho, que a associação com a chuva seja negativa.
# No modelo sem interação, gamma = 0 por definição, não por estimação.
# Os coeficientes descrevem associações preditivas, não efeitos causais.


# As previsões fora de [0,1] ocorreram apenas em GWETROOT de Belém e Manaus.
# Foram 120 casos sem interação e 87 com interação, todos acima de 1.
# A melhora das métricas não eliminou essa limitação dos modelos.


# Gráficos de comparação e previsão

dir.create(pasta_resultados, recursive = TRUE, showWarnings = FALSE)
dir.create(pasta_figuras, recursive = TRUE, showWarnings = FALSE)

ganhos <- selecionados |>
  select(localidade, variavel, familia, RMSE = ganho_RMSE_AR, MAE = ganho_MAE_AR) |>
  pivot_longer(c(RMSE, MAE), names_to = "metrica", values_to = "ganho") |>
  mutate(metrica = factor(metrica, levels = c("RMSE", "MAE")))

grafico_metricas <- ggplot(ganhos, aes(ganho, localidade, color = familia)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(position = position_dodge(width = 0.5), size = 2.2) +
  facet_wrap(vars(metrica, variavel), ncol = 2, scales = "free_x") +
  scale_color_manual(values = cores) +
  labs(x = "Redução do erro em relação ao AR (%)", y = NULL, color = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

curvas <- bind_rows(previsoes |> semi_join(selecionados, by = chaves_modelo),
                    previsoes_ar |> mutate(familia = "AR")) |>
  filter(data_previsao >= inicio_grafico, data_previsao <= fim_grafico,
         localidade %in% localidades_grafico)
observados <- curvas |> distinct(localidade, variavel, data_previsao, observado)

grafico_previsao <- ggplot() +
  geom_line(data = curvas, aes(data_previsao, previsto, color = familia), linewidth = 0.5) +
  geom_line(data = observados, aes(data_previsao, observado, color = "Observado"),
            linewidth = 0.4, linetype = "dashed") +
  facet_wrap(vars(localidade, variavel), ncol = 2, scales = "free_y") +
  scale_x_date(date_breaks = "1 month", date_labels = "%m/%Y", expand = expansion(mult = 0)) +
  scale_color_manual(values = cores, breaks = c("Observado", "AR", "AR + P", "AR + P + interação")) +
  labs(x = NULL, y = "Umidade do solo", color = NULL,
       subtitle = paste("Previsões de um dia à frente:", format(inicio_grafico, "%d/%m/%Y"),
                        "a", format(fim_grafico, "%d/%m/%Y"))) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

# O recorte de outubro a dezembro facilita enxergar as diferenças entre as curvas.
# As métricas usam os três anos completos; o recorte altera apenas o gráfico.

if (interactive()) {
  print(grafico_metricas)
  print(grafico_previsao)
}


# Exportação para conferência e redação do artigo

write_csv(previsoes, file.path(pasta_resultados, "previsoes_validacao.csv"))
write_csv(metricas, file.path(pasta_resultados, "metricas_validacao.csv"))
write_csv(metricas_anuais, file.path(pasta_resultados, "metricas_validacao_anuais.csv"))
write_csv(selecionados, file.path(pasta_resultados, "selecionados_validacao.csv"))
write_csv(coeficientes_precipitacao, file.path(pasta_resultados, "coeficientes_precipitacao.csv"))
write_csv(falhas, file.path(pasta_resultados, "falhas.csv"))
saveRDS(ajustes_selecionados, file.path(pasta_resultados, "ajustes_selecionados.rds"))

ggsave(file.path(pasta_figuras, "precipitacao_validacao.pdf"),
       grafico_metricas, width = 9, height = 6.5)
ggsave(file.path(pasta_figuras, "precipitacao_observado_previsto.pdf"),
       grafico_previsao, width = 10, height = 8)

# A comparação avalia a chuva mantendo a estrutura temporal escolhida no AR.
# Os resultados de 2025 não entram nos ajustes nem na seleção desta etapa.
