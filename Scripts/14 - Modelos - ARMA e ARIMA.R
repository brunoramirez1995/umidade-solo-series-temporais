
# 14 - Regressão ARMA/ARIMA

# Pacotes e resultados dos scripts 12 e 13

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

source("Scripts/apoio_modelagem.R", encoding = "UTF-8")


arquivo_dados <- "Dados/nasa_power_diario.csv"
pasta_resultados <- "resultados/modelagem/arma_arima"

chaves_serie <- c("localidade_id", "variavel")
chaves_modelo <- c(chaves_serie, "modelo")

referencia <- read_csv(
  "resultados/modelagem/referencias/metricas_validacao.csv",
  show_col_types = FALSE
) |>
  filter(modelo == "persistencia") |>
  select(all_of(chaves_serie), RMSE_persistencia = RMSE, MAE_persistencia = MAE)



referencia_ar <- read_csv(
  "resultados/modelagem/autorregressivos/selecionados_validacao.csv",
  col_types = cols(K = col_integer(), tendencia = col_logical()),
  show_col_types = FALSE
) |>
  select(all_of(chaves_serie), K, tendencia, modelo_AR = modelo,
         RMSE_AR = RMSE, MAE_AR = MAE)


# Base diária e calendário


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
  mutate(tempo = as.numeric(data - inicio_treinamento) / 365.25)



# Cada série precisa de uma configuração de calendário do script 13

for (k in seq_len(max(referencia_ar$K))) {
  dados[[paste0("sen_", k)]] <- sin(2 * pi * k * dados$tempo)
  dados[[paste0("cos_", k)]] <- cos(2 * pi * k * dados$tempo)
}


# Sete candidatos por série
# K e tendência vêm do AR; todos os coeficientes serão estimados novamente


ordens <- tibble::tribble(
  ~p, ~d, ~q,
   1,  0,  1,
   2,  0,  1,
   1,  0,  2,
   0,  1,  1,
   1,  1,  1,
   2,  1,  1,
   1,  1,  2
)

configuracoes <- referencia_ar |>
  select(all_of(chaves_serie), K, tendencia) |>
  crossing(ordens) |>
  mutate(
    modelo = paste0("ARIMA_", p, "_", d, "_", q, "_K", K,
                    if_else(tendencia, "_T", "")),
    n_parametros = p + q + 2L * K + as.integer(tendencia) +
      as.integer(d == 0) + 1L
  )

# A contagem de parâmetros inclui a variância das inovações
# Os termos AR e MA pertencem ao erro da regressão; não são colunas em xreg.


# Previsões diárias com coeficientes fixos durante cada ano

prever_diariamente <- function(ajuste, observado, xreg) {
  
  parte_regressao <- as.numeric(xreg %*% coef(ajuste)[colnames(xreg)])
  
  if ("intercept" %in% names(coef(ajuste))) {
    parte_regressao <- parte_regressao + coef(ajuste)["intercept"]
  }

  estado <- ajuste$model
  previsto <- numeric(length(observado))
  
  for (t in seq_along(observado)) {
    previsto[t] <- parte_regressao[t] +
      as.numeric(stats::KalmanForecast(1L, estado)$pred)
    atualizacao <- stats::KalmanRun(
      observado[t] - parte_regressao[t], estado, nit = -1L, update = TRUE
    )
    estado <- attr(atualizacao, "mod")
  }
  previsto
}

# Primeiro se prevê; depois se incorpora a observação ao estado do ARIMA.
# nit = -1 permite avançar a covariância do filtro a cada novo dia.


# Ajustes anuais e previsões na validação

series <- dados |> group_by(localidade_id, variavel) |> group_split()

lista_previsoes <- list()
lista_ajustes <- list()

falhas <- tibble::tibble(
  localidade_id = character(), variavel = character(), modelo = character(),
  ano = integer(), motivo = character()
)

# O estado do ARIMA resume o histórico necessário para calcular as previsões.
# Em cada dia, primeiro se prevê e depois se incorpora o valor observado.
# Os coeficientes permanecem fixos até a reestimação no ano seguinte.


# São 12 séries, sete candidatos e três anos de validação.
# Cada ajuste utiliza somente os dados disponíveis até o ano anterior.

for (serie in series) {
  
  message(serie$localidade[1], " - ", serie$variavel[1])
  candidatas <- configuracoes |>
    filter(localidade_id == serie$localidade_id[1], variavel == serie$variavel[1])

  for (i in seq_len(nrow(candidatas))) {
    config <- candidatas[i, ]
    entradas <- c(
      paste0("sen_", seq_len(config$K)),
      paste0("cos_", seq_len(config$K)),
      if (config$tendencia) "tempo"
    )

    for (j in seq_len(nrow(calendario_validacao))) {
      janela <- calendario_validacao[j, ]
      treino <- serie |> filter(data <= janela$fim_ajuste)
      avaliacao <- serie |>
        filter(data >= janela$inicio_avaliacao, data <= janela$fim_avaliacao)

      resultado <- tryCatch({
        x_treino <- as.matrix(select(treino, all_of(entradas)))
        x_avaliacao <- as.matrix(select(avaliacao, all_of(entradas)))

        ajuste <- stats::arima(
          treino$observado, order = c(config$p, config$d, config$q),
          seasonal = c(0, 0, 0), xreg = x_treino,
          include.mean = config$d == 0, method = "ML",
          SSinit = "Rossignol2011", optim.control = list(maxit = 1000)
        )
        if (ajuste$code != 0L) stop("Ajuste sem convergência: código ", ajuste$code, ".")
        if (any(!is.finite(coef(ajuste)))) stop("Coeficientes não finitos.")
        
        # Verifica a matriz de variância dos coeficientes
        
        if (any(!is.finite(ajuste$var.coef)) ||
            inherits(try(chol(ajuste$var.coef), silent = TRUE), "try-error")) {
          stop("Covariância dos coeficientes não positiva definida.")
        }
        

        previsto <- prever_diariamente(ajuste, avaliacao$observado, x_avaliacao)
        if (any(!is.finite(previsto))) stop("Previsões não finitas.")

        list(
          ajuste = ajuste,
          previsoes = avaliacao |>
            transmute(
              localidade_id, localidade, clima, variavel,
              data_origem = data - 1, data_previsao = data,
              ano = janela$ano, modelo = config$modelo,
              observado, previsto = .env$previsto, erro = observado - previsto
            )
        )
      }, error = function(e) e)

      if (inherits(resultado, "error")) {
        falhas <- bind_rows(falhas, tibble::tibble(
          localidade_id = serie$localidade_id[1], variavel = serie$variavel[1],
          modelo = config$modelo, ano = janela$ano, motivo = conditionMessage(resultado)
        ))
        message("Falha: ", config$modelo, " / ", janela$ano)
      } else {
        lista_previsoes[[length(lista_previsoes) + 1L]] <- resultado$previsoes
        lista_ajustes[[length(lista_ajustes) + 1L]] <- tibble::tibble(
          localidade_id = serie$localidade_id[1], variavel = serie$variavel[1],
          modelo = config$modelo, ano = janela$ano,
          fim_ajuste = janela$fim_ajuste,
          dados_treino = list(select(treino, data, observado)),
          ajuste = list(resultado$ajuste)
        )
      }
    }
  }
}

# Em d = 1, arima() diferencia a relação completa, inclusive os regressores
# A tendência linear representa a deriva; não há outro intercepto em níveis

dir.create(pasta_resultados, recursive = TRUE, showWarnings = FALSE)

write_csv(falhas, file.path(pasta_resultados, "falhas.csv"))



previsoes <- bind_rows(lista_previsoes)
chaves <- c("localidade_id", "localidade", "clima", "variavel", "modelo")


# Métricas e seleção de p, d e q

metricas <- previsoes |>
  group_by(across(all_of(chaves))) |>
  summarise(
    calcular_metricas(observado, previsto),
    n_fora = sum(previsto < 0 | previsto > 1),
    excesso_maximo = max(pmax(-previsto, previsto - 1, 0)),
    .groups = "drop"
  ) |>
  left_join(configuracoes, by = chaves_modelo) |>
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
  ungroup() |>
  left_join(referencia, by = chaves_serie) |>
  left_join(select(referencia_ar, -K, -tendencia), by = chaves_serie) |>
  mutate(
    ganho_RMSE = 100 * (1 - RMSE / RMSE_persistencia),
    ganho_MAE = 100 * (1 - MAE / MAE_persistencia),
    ganho_RMSE_AR = 100 * (1 - RMSE / RMSE_AR),
    ganho_MAE_AR = 100 * (1 - MAE / MAE_AR)
  )

# O RMSE usa todos os dias de 2022 a 2024. Candidatos incompletos não entram
# As previsões não são truncadas; 2025 permanece fora desta etapa


previsoes_selecionadas <- previsoes |>
  semi_join(selecionados, by = chaves_modelo)

ajustes_selecionados <- bind_rows(lista_ajustes) |>
  semi_join(selecionados, by = chaves_modelo) |>
  left_join(configuracoes, by = chaves_modelo)


# Conferência e exportação



message("Séries com configuração completa: ", nrow(selecionados), " de ", length(series))

print(selecionados |>
        select(localidade, variavel, modelo, p, d, q, K, tendencia,
               RMSE, MAE, vies, ganho_RMSE, ganho_MAE,
               ganho_RMSE_AR, ganho_MAE_AR, n_fora, excesso_maximo),
      n = Inf, width = Inf)


# Nesta rodada, o RMSE foi menor que o da persistência nas 12 séries,
# mas o MAE aumentou em todas elas.

# Nos dias sem mudança nos valores registrados, a persistência acerta
# exatamente. Os modelos introduzem pequenos desvios nesses dias.
# Nos dias com mudança, os modelos reduziram os erros absolutos e quadráticos.

# Frente ao AR, houve redução do RMSE em sete séries, sempre inferior a 1%.
# Os resultados não mostram uma vantagem ampla dos ARMA/ARIMA.

# Em Brasília/GWETTOP, foi selecionado o ARMA(1,2), com K = 1 e tendência.
# O ARIMA(2,1,1) foi excluído pelo problema de covariância no ajuste de 2024.



print(metricas_anuais |> semi_join(selecionados, by = chaves_modelo),
      n = Inf, width = Inf)

# O RMSE foi menor que o da persistência em 35 das 36 combinações
# entre série e ano.
# A exceção foi GWETROOT de Petrolina em 2023, com aumento de 1,64%.
# Portanto, o ganho observado no período completo não ocorreu em todos os anos.



print(falhas, n = Inf, width = Inf)

# Dos 252 ajustes anuais, 15 foram rejeitados pelo critério de covariância,
# envolvendo oito configurações.
# Restaram 76 configurações completas, avaliadas nos 1.096 dias de validação.
# As configurações parciais ficaram fora da seleção.


write_csv(previsoes, file.path(pasta_resultados, "previsoes_validacao.csv"))
write_csv(metricas, file.path(pasta_resultados, "metricas_validacao.csv"))
write_csv(metricas_anuais, file.path(pasta_resultados, "metricas_validacao_anuais.csv"))
write_csv(selecionados, file.path(pasta_resultados, "selecionados_validacao.csv"))
saveRDS(ajustes_selecionados, file.path(pasta_resultados, "ajustes_selecionados.rds"))

# O RDS guarda os três ajustes anuais de cada vencedor e seus dados de treino

# Ganhos positivos indicam redução do erro frente à referência indicada

# A comparação com o AR não isola o efeito de MA: formulação e estimação mudam


# Comparação dos erros e curvas de previsão



# Cada ponto previsto usa apenas as informações disponíveis no dia anterior.
# O gráfico reúne previsões sucessivas de um dia à frente.
# A comparação numérica com as referências é necessária porque
# a sobreposição das curvas pode esconder diferenças pequenas.


if (nrow(selecionados) > 0) {
  dir.create("figuras/principais", recursive = TRUE, showWarnings = FALSE)
  comparacao <- selecionados |>
    select(localidade, variavel, ganho_RMSE, ganho_RMSE_AR) |>
    pivot_longer(starts_with("ganho_"), names_to = "referencia", values_to = "ganho") |>
    mutate(referencia = recode(
      referencia, ganho_RMSE = "Persistência", ganho_RMSE_AR = "AR com harmônicos"
    ))

  grafico <- ggplot(comparacao, aes(x = ganho, y = localidade)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.5, color = "#2166AC") +
    facet_grid(referencia ~ variavel) +
    labs(x = "Redução do RMSE em relação à referência (%)", y = NULL) +
    theme_minimal(base_size = 12)

  ggsave("figuras/principais/arma_arima_validacao.pdf",
         grafico, width = 9, height = 6.5)
  if (interactive()) print(grafico)

  curvas <- previsoes_selecionadas |>
    filter(ano == 2024, localidade_id %in% c("belem_pa", "brasilia_df", "petrolina_pe")) |>
    select(data_previsao, localidade, variavel, observado, previsto) |>
    pivot_longer(c(observado, previsto), names_to = "curva", values_to = "umidade")

  if (nrow(curvas) > 0) {
    grafico_previsao <- ggplot(curvas, aes(data_previsao, umidade, color = curva)) +
      geom_line(linewidth = 0.4) +
      facet_wrap(vars(localidade, variavel), ncol = 2, scales = "free_y") +
      scale_x_date(date_breaks = "3 months", date_labels = "%m/%Y",
                   expand = expansion(mult = 0)) +
      scale_color_manual(values = c(observado = "grey25", previsto = "#2166AC"),
                         labels = c(observado = "Observado", previsto = "Previsto")) +
      labs(x = NULL, y = "Umidade do solo", color = NULL,
           subtitle = "Previsões de um dia à frente em 2024") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")

    ggsave("figuras/principais/arma_arima_observado_previsto.pdf",
           grafico_previsao, width = 10, height = 8)
    if (interactive()) print(grafico_previsao)
  }
}
