# Modelos AR com harmônicos: previsão de um dia à frente

# Pacotes e configurações

library(tidyverse)
source("Scripts/apoio_modelagem.R", encoding = "UTF-8")

arquivo_dados <- "Dados/nasa_power_diario.csv"
pasta_resultados <- "resultados/modelagem/autorregressivos"


# Leitura e construção das entradas

dados <- read_csv(
  arquivo_dados,
  col_types = cols_only(
    localidade_id = col_character(), localidade = col_character(),
    clima = col_character(), data = col_date(),
    GWETTOP = col_double(), GWETROOT = col_double()
  ),
  show_col_types = FALSE
) |>
  
  filter(data >= inicio_treinamento, data <= fim_validacao) |>
  
  pivot_longer(all_of(alvos), names_to = "variavel", values_to = "observado") |>
  
  arrange(localidade_id, variavel, data) |>
  
  group_by(localidade_id, variavel) |>
  mutate(
    lag_1 = dplyr::lag(observado, 1),
    lag_2 = dplyr::lag(observado, 2),
    lag_3 = dplyr::lag(observado, 3)
  ) |>
  
  ungroup() |>
  mutate(tempo = as.numeric(data - inicio_treinamento) / 365.25)


summary(dados)


# O calendário é contínuo, inclusive em anos bissextos.
# Os harmônicos correspondem à data prevista; os lags, aos dias anteriores.

for (k in 1:3) {
  dados[[paste0("sen_", k)]] <- sin(2 * pi * k * dados$tempo)
  dados[[paste0("cos_", k)]] <- cos(2 * pi * k * dados$tempo)
}


# Candidatos definidos na metodologia: 18 configurações por série

configuracoes <- expand_grid(p = 1:3, K = 1:3, tendencia = c(FALSE, TRUE)) |>
  mutate(
    modelo = paste0("AR", p, "_K", K, if_else(tendencia, "_T", "")),
    n_parametros = 2L + p + 2L * K + as.integer(tendencia)
  )

# A contagem inclui os coeficientes da regressão e a variância dos erros.
# O sufixo _T identifica a inclusão de tendência linear.


# Ajustes anuais e previsões diárias

series <- dados |>
  group_by(localidade_id, variavel) |>
  group_split()

lista_previsoes <- list()

falhas <- tibble(
  localidade_id = character(), variavel = character(), modelo = character(),
  ano = integer(), motivo = character()
)

for (serie in series) {
  message(serie$localidade[1], " - ", serie$variavel[1])

  for (i in seq_len(nrow(configuracoes))) {
    config <- configuracoes[i, ]
    entradas <- c(
      paste0("lag_", seq_len(config$p)),
      paste0("sen_", seq_len(config$K)),
      paste0("cos_", seq_len(config$K)),
      if (config$tendencia) "tempo"
    )
    formula_ar <- reformulate(entradas, response = "observado")

    for (j in seq_len(nrow(calendario_validacao))) {
      janela <- calendario_validacao[j, ]
      treino <- serie |> filter(data <= janela$fim_ajuste)
      avaliacao <- serie |>
        filter(data >= janela$inicio_avaliacao, data <= janela$fim_avaliacao)

      resultado <- tryCatch({
        ajuste <- lm(formula_ar, data = treino, na.action = na.omit)
        if (anyNA(coef(ajuste))) stop("Coeficientes não identificados.")

        previsto <- as.numeric(predict(ajuste, newdata = avaliacao))
        if (any(!is.finite(previsto))) stop("Previsões não finitas.")

        avaliacao |>
          transmute(
            localidade_id, localidade, clima, variavel,
            data_origem = data - 1, data_previsao = data,
            ano = janela$ano, modelo = config$modelo,
            observado, previsto = .env$previsto, erro = observado - previsto
          )
      }, error = function(e) e)

      # Uma falha de ajuste é registrada, sem interromper os outros candidatos.
      if (inherits(resultado, "error")) {
        falhas <- bind_rows(falhas, tibble(
          localidade_id = serie$localidade_id[1], variavel = serie$variavel[1],
          modelo = config$modelo, ano = janela$ano,
          motivo = conditionMessage(resultado)
        ))
        message("Falha: ", config$modelo, " / ", janela$ano)
      } else {
        lista_previsoes[[length(lista_previsoes) + 1L]] <- resultado
      }
    }
  }
}

# Os parâmetros ficam fixos durante cada ano. As entradas são atualizadas
# diariamente com a umidade observada até o dia anterior à previsão.
# As primeiras p observações do treinamento não têm todos os lags necessários.

previsoes <- bind_rows(lista_previsoes)
chaves <- c("localidade_id", "localidade", "clima", "variavel", "modelo")


# Métricas e seleção preliminar dentro da família AR

metricas <- previsoes |>
  group_by(across(all_of(chaves))) |>
  summarise(
    calcular_metricas(observado, previsto),
    n_fora = sum(previsto < 0 | previsto > 1),
    excesso_maximo = max(pmax(-previsto, previsto - 1, 0)),
    .groups = "drop"
  ) |>
  left_join(configuracoes, by = "modelo") |>
  mutate(
    completo = n == as.integer(fim_validacao - inicio_validacao) + 1L,
    percentual_fora = 100 * n_fora / n
  )

metricas_anuais <- previsoes |>
  group_by(across(all_of(c(chaves, "ano")))) |>
  summarise(calcular_metricas(observado, previsto), .groups = "drop")

selecionados <- metricas |>
  filter(completo) |>
  group_by(localidade_id, variavel) |>
  arrange(RMSE, n_parametros, modelo, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup()

# Somente configurações com todo o período previsto entram na seleção.
# O menor RMSE define o candidato; em empate, vence o menos parametrizado.
# As previsões são avaliadas como saem do modelo, sem limitar a [0,1].


# Salvar os ajustes selecionados para os diagnósticos

# São refeitos somente os 36 ajustes selecionados para salvamento.
# Isso evita manter os 648 ajustes da grade na memória.

lista_ajustes <- list()

for (i in seq_len(nrow(selecionados))) {
  config <- selecionados[i, ]
  entradas <- c(
    paste0("lag_", seq_len(config$p)),
    paste0("sen_", seq_len(config$K)),
    paste0("cos_", seq_len(config$K)),
    if (config$tendencia) "tempo"
  )
  formula_ar <- reformulate(entradas, response = "observado")

  for (j in seq_len(nrow(calendario_validacao))) {
    janela <- calendario_validacao[j, ]
    treino <- dados |>
      filter(localidade_id == config$localidade_id, variavel == config$variavel,
             data >= janela$inicio_ajuste, data <= janela$fim_ajuste) |>
      drop_na(observado, all_of(entradas))

    ajuste <- lm(formula_ar, data = treino, na.action = na.omit)

    lista_ajustes[[length(lista_ajustes) + 1L]] <- tibble(
      localidade_id = config$localidade_id, variavel = config$variavel,
      modelo = config$modelo, p = config$p, K = config$K,
      tendencia = config$tendencia, ano = janela$ano,
      inicio_ajuste = janela$inicio_ajuste, fim_ajuste = janela$fim_ajuste,
      dados_treino = list(select(treino, data, observado)),
      ajuste = list(ajuste)
    )
  }
}

ajustes_selecionados <- bind_rows(lista_ajustes)


# Comparação com a persistência calculada no script 12

referencia <- read_csv(
  "resultados/modelagem/referencias/metricas_validacao.csv",
  show_col_types = FALSE
) |>
  filter(modelo == "persistencia") |>
  select(localidade_id, variavel,
         RMSE_persistencia = RMSE, MAE_persistencia = MAE)

selecionados <- selecionados |>
  left_join(referencia, by = c("localidade_id", "variavel")) |>
  mutate(
    ganho_RMSE = 100 * (1 - RMSE / RMSE_persistencia),
    ganho_MAE = 100 * (1 - MAE / MAE_persistencia)
  )

message("Séries com configuração completa: ", nrow(selecionados), " de ", length(series))

print(selecionados |>
        select(localidade, variavel, modelo, p, K, tendencia,
               RMSE, MAE, vies, ganho_RMSE, ganho_MAE,
               n_fora, excesso_maximo), n = Inf, width = Inf)

print(metricas_anuais |>
        semi_join(selecionados, by = c("localidade_id", "variavel", "modelo")),
      n = Inf, width = Inf)
print(falhas, n = Inf, width = Inf)

# Ganho positivo indica redução do erro em relação à persistência.
# Esta seleção ainda será examinada pelos diagnósticos dos resíduos.
# K começa em 1: este bloco compara AR com harmônicos, sem isolar seu efeito.


# Comentários sobre os resultados

# Os AR selecionados reduziram o RMSE nas 12 séries frente à persistência.
# Os ganhos variaram de 1,98% a 11,42%.

# O MAE aumentou em 11 séries. A única redução foi praticamente nula:
# 0,05% em GWETTOP de Petrolina.

# Nos dias com mudança registrada na umidade, os AR tiveram menor MAE.
# Nos dias sem mudança, a persistência acertava exatamente, enquanto
# os AR introduziam pequenos erros. Isso ajuda a explicar a diferença
# entre os resultados de RMSE e MAE.

# O RMSE melhorou em 34 das 36 comparações anuais.
# As exceções foram Petrolina/GWETROOT em 2023 e Uberlândia/GWETTOP em 2024.

# Os candidatos mais bem classificados tiveram resultados muito próximos.
# A seleção não demonstra que uma ordem específica seja indispensável.

# Houve previsões acima de 1 em GWETROOT de Belém e Manaus.
# As previsões foram mantidas sem truncamento, conforme o protocolo.

# Estes resultados orientam a seleção. A avaliação final ainda depende
# dos diagnósticos dos resíduos e do desempenho no teste de 2025.




# Gráfico e exportação

grafico <- ggplot(selecionados, aes(x = ganho_RMSE, y = localidade)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5, color = "#2166AC") +
  facet_wrap(~ variavel) +
  labs(x = "Redução do RMSE em relação à persistência (%)", y = NULL) +
  theme_minimal(base_size = 12)

if (interactive()) print(grafico)

dir.create(pasta_resultados, recursive = TRUE, showWarnings = FALSE)
dir.create("figuras/principais", recursive = TRUE, showWarnings = FALSE)

write_csv(previsoes, file.path(pasta_resultados, "previsoes_validacao.csv"))
write_csv(metricas, file.path(pasta_resultados, "metricas_validacao.csv"))
write_csv(metricas_anuais, file.path(pasta_resultados, "metricas_validacao_anuais.csv"))
write_csv(selecionados, file.path(pasta_resultados, "selecionados_validacao.csv"))
write_csv(falhas, file.path(pasta_resultados, "falhas.csv"))
saveRDS(ajustes_selecionados, file.path(pasta_resultados, "ajustes_selecionados.rds"))

ggsave("figuras/principais/ar_harmonicos_validacao.pdf",
       grafico, width = 9, height = 4.5)
