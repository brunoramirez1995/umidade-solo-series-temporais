# 18 - Comparações e diagnósticos dos modelos


# Pacotes e definições dos scripts anteriores

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

source("Scripts/apoio_modelagem.R", encoding = "UTF-8")

pasta_resultados <- "resultados/modelagem/comparacoes"
pasta_figuras <- "figuras/principais"
pasta_diagnosticos <- "figuras/diagnosticos"

dir.create(pasta_resultados, recursive = TRUE, showWarnings = FALSE)
dir.create(pasta_figuras, recursive = TRUE, showWarnings = FALSE)
dir.create(pasta_diagnosticos, recursive = TRUE, showWarnings = FALSE)

ano_grafico <- 2024
localidades_grafico <- c("Belém", "Brasília", "Petrolina")
familias <- c("Persistência", "AR", "ARMA/ARIMA", "AR + P", "AR + P + interação",
              "SETAR", "LSTAR")
cores <- c("Persistência" = "#D55E00", "AR" = "#009E73",
           "ARMA/ARIMA" = "#2166AC", "AR + P" = "#CC79A7",
           "AR + P + interação" = "#E69F00", "SETAR" = "#7B3294",
           "LSTAR" = "#A65628", "Observado" = "grey25")

# As figuras lineares são mantidas; as gerais incluem todas as famílias.
comparacoes <- list(lineares = familias[1:3], modelos = familias)

chaves_modelo <- c("localidade_id", "variavel", "modelo")
chaves_serie <- c("localidade_id", "localidade", "clima", "variavel")


# Leitura das previsões e das configurações selecionadas

selecionados_ar <- read_csv(
  "resultados/modelagem/autorregressivos/selecionados_validacao.csv",
  show_col_types = FALSE
)

selecionados_arma <- read_csv(
  "resultados/modelagem/arma_arima/selecionados_validacao.csv",
  show_col_types = FALSE
)

selecionados_precipitacao <- read_csv(
  "resultados/modelagem/precipitacao/selecionados_validacao.csv",
  show_col_types = FALSE
)

selecionados_setar <- read_csv(
  "resultados/modelagem/setar/selecionados_validacao.csv",
  show_col_types = FALSE
)

selecionados_lstar <- read_csv(
  "resultados/modelagem/lstar/selecionados_validacao.csv",
  show_col_types = FALSE
)

referencias <- read_csv(
  "resultados/modelagem/referencias/previsoes_validacao.csv",
  show_col_types = FALSE
) |>
  filter(modelo == "persistencia") |>
  mutate(familia = "Persistência")

previsoes_ar <- read_csv(
  "resultados/modelagem/autorregressivos/previsoes_validacao.csv",
  show_col_types = FALSE
) |>
  semi_join(selecionados_ar, by = chaves_modelo) |>
  mutate(familia = "AR")

previsoes_arma <- read_csv(
  "resultados/modelagem/arma_arima/previsoes_validacao.csv",
  show_col_types = FALSE
) |>
  semi_join(selecionados_arma, by = chaves_modelo) |>
  mutate(familia = "ARMA/ARIMA")

previsoes_precipitacao <- read_csv(
  "resultados/modelagem/precipitacao/previsoes_validacao.csv",
  show_col_types = FALSE
) |>
  semi_join(selecionados_precipitacao, by = chaves_modelo)

previsoes_setar <- read_csv(
  "resultados/modelagem/setar/previsoes_validacao.csv",
  show_col_types = FALSE
) |>
  semi_join(selecionados_setar, by = chaves_modelo) |>
  mutate(familia = "SETAR")

previsoes_lstar <- read_csv(
  "resultados/modelagem/lstar/previsoes_validacao.csv",
  show_col_types = FALSE
) |>
  semi_join(selecionados_lstar, by = chaves_modelo) |>
  mutate(familia = "LSTAR")

previsoes <- bind_rows(referencias, previsoes_ar, previsoes_arma,
                       previsoes_precipitacao, previsoes_setar, previsoes_lstar) |>
  mutate(familia = factor(familia, levels = familias)) |>
  arrange(localidade_id, variavel, familia, data_previsao)

# A comparação exige as mesmas datas e observações nas sete famílias.

conferencia <- previsoes |>
  group_by(localidade_id, variavel, data_previsao) |>
  summarise(n = n(), familias = n_distinct(familia),
            observados = n_distinct(observado), .groups = "drop")

if (anyNA(previsoes$familia) ||
    any(conferencia$n != length(familias) |
        conferencia$familias != length(familias) | conferencia$observados != 1)) {
  stop("As previsões das famílias não correspondem aos mesmos dados.")
}


# Métricas no período completo e em cada ano

resumir_previsoes <- function(base, grupos) {

  metricas <- base |>
    group_by(across(all_of(c(grupos, "familia", "modelo")))) |>
    summarise(
      calcular_metricas(observado, previsto),
      n_fora = sum(previsto < 0 | previsto > 1),
      percentual_fora = 100 * mean(previsto < 0 | previsto > 1),
      excesso_maximo = max(pmax(-previsto, previsto - 1, 0)),
      .groups = "drop"
    )

  persistencia <- metricas |>
    filter(familia == "Persistência") |>
    select(all_of(grupos), RMSE_persistencia = RMSE, MAE_persistencia = MAE)

  referencia_ar <- metricas |>
    filter(familia == "AR") |>
    select(all_of(grupos), RMSE_AR = RMSE, MAE_AR = MAE)

  metricas |>
    left_join(persistencia, by = grupos) |>
    left_join(referencia_ar, by = grupos) |>
    mutate(
      ganho_RMSE = 100 * (1 - RMSE / RMSE_persistencia),
      ganho_MAE = 100 * (1 - MAE / MAE_persistencia),
      ganho_RMSE_AR = 100 * (1 - RMSE / RMSE_AR),
      ganho_MAE_AR = 100 * (1 - MAE / MAE_AR)
    )
}

metricas <- resumir_previsoes(previsoes, chaves_serie)
metricas_anuais <- resumir_previsoes(previsoes, c(chaves_serie, "ano"))

print(metricas |>
        select(localidade, variavel, familia, modelo, n, RMSE, MAE,
               ganho_RMSE, ganho_MAE, ganho_RMSE_AR, ganho_MAE_AR,
               vies, n_fora, excesso_maximo),
      n = Inf, width = Inf)

print(metricas_anuais |>
        select(localidade, variavel, familia, ano, RMSE, MAE,
               ganho_RMSE, ganho_MAE, ganho_RMSE_AR, ganho_MAE_AR, vies, n_fora),
      n = Inf, width = Inf)

# ganho_RMSE e ganho_MAE usam a persistência; as colunas com _AR usam o AR.
# Ganhos positivos indicam redução do erro frente à referência correspondente.
# O viés é observado - previsto: positivo indica subestimação média.
# As previsões são avaliadas como foram exportadas, sem limitar a [0,1].

# Resultados de 2022 a 2024: os ARMA/ARIMA reduziram o RMSE da persistência
# nas 12 séries. Frente ao AR, houve redução em sete, sempre inferior a 1%
# no período completo. Nas comparações anuais, venceram em 25 de 36 casos.

# Entre persistência, AR e ARMA/ARIMA, a persistência teve o menor MAE em 11 séries.
# Em GWETTOP de Petrolina, o AR teve uma vantagem de apenas 0,05%.

# Houve previsões acima de 1 apenas em GWETROOT de Belém e Manaus.
# Em Manaus, todas ocorreram em 2022: 45 no AR e 54 no ARMA/ARIMA.
# Nesse ano, correspondem a 12,33% e 14,79% dos dias, respectivamente.

# Em Uberlândia/GWETTOP de 2024, os dois modelos superestimaram a umidade, em média.
# O viés foi de -0,00501 no AR e -0,00456 no ARMA/ARIMA.


# A chuva reduziu o RMSE de GWETROOT nas seis localidades.
# Em GWETTOP, o ganho apareceu principalmente com a interação, sobretudo em Petrolina.
# A interação melhorou o RMSE frente ao AR + P em 11 séries e o MAE em oito.
# As previsões acima de 1 passaram de 120 no AR + P para 87 com interação.




# Entre as sete famílias, a interação teve o menor RMSE em oito séries.
# O LSTAR venceu em GWETTOP de Brasília e Manaus; o SETAR, em Uberlândia.
# Em Uberlândia/GWETROOT, o menor RMSE continuou sendo do AR + P.
# Pelo MAE, a persistência continuou sendo a melhor em três séries.
# A interação venceu pelo RMSE em sete séries em 2022, cinco em 2023
# e 11 em 2024. O LSTAR venceu em três, duas e uma, respectivamente.
# A comparação principal considera os três anos.
# Houve menos previsões acima de 1 com interação, mas o maior excesso
# ainda foi superior ao observado no AR e no ARMA/ARIMA.



# O SETAR reduziu o RMSE do AR nas 12 séries e o MAE em oito.
# Os maiores ganhos de RMSE ocorreram em GWETTOP de Brasília (7,06%)
# e Petrolina (8,01%). Em GWETROOT, todos os ganhos ficaram abaixo de 0,8%.
# Foram obtidas 50 previsões acima de 1 no SETAR: quatro em Belém e 46 em Manaus,
# todas em GWETROOT. O maior valor previsto foi 1,00555.


# O LSTAR reduziu o RMSE do AR nas 12 séries e o do SETAR em oito.
# Os ganhos adicionais frente ao SETAR foram mais claros em GWETTOP.
# Em GWETROOT, as diferenças de RMSE ficaram abaixo de 0,4% em módulo.
# A redução do RMSE não foi acompanhada por redução do MAE em todos os casos.
# Houve 29 previsões acima de 1: dez em Belém e 19 em Manaus, em GWETROOT.
# O maior valor previsto foi 1,00247. Os limites da transição são avaliados no script 17.


# Comparação gráfica de RMSE e MAE

for (comparacao in names(comparacoes)) {
  familias_grafico <- comparacoes[[comparacao]]
  linhas_legenda <- if (comparacao == "lineares") 1 else 2

  metricas_longas <- metricas |>
    filter(familia %in% familias_grafico) |>
    pivot_longer(c(RMSE, MAE), names_to = "metrica", values_to = "valor") |>
    mutate(metrica = factor(metrica, levels = c("RMSE", "MAE")))

  grafico_metricas <- ggplot(metricas_longas,
                            aes(valor, localidade, color = familia)) +
    geom_point(position = position_dodge(width = 0.55), size = 2.2) +
    facet_wrap(vars(metrica, variavel), ncol = 2, scales = "free_x") +
    scale_color_manual(values = cores) +
    labs(x = "Erro na validação de 2022 a 2024", y = NULL, color = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = linhas_legenda, byrow = TRUE))

  ggsave(file.path(pasta_figuras, paste0("comparacao_", comparacao, "_validacao.pdf")),
         grafico_metricas, width = 9, height = 6.5)
  if (interactive()) print(grafico_metricas)

  metricas_anuais_grafico <- metricas_anuais |>
    filter(familia %in% familias_grafico)

  for (metrica in c("RMSE", "MAE")) {
    grafico_anual <- ggplot(metricas_anuais_grafico,
                           aes(ano, .data[[metrica]], color = familia)) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 1.5) +
      facet_wrap(vars(localidade, variavel), ncol = 2, scales = "free_y") +
      scale_x_continuous(breaks = 2022:2024) +
      scale_color_manual(values = cores) +
      labs(x = NULL, y = metrica, color = NULL) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom") +
      guides(color = guide_legend(nrow = linhas_legenda, byrow = TRUE))

    ggsave(file.path(pasta_figuras, paste0("comparacao_", comparacao, "_", metrica, "_anual.pdf")),
           grafico_anual, width = 9, height = 11)
    if (interactive()) print(grafico_anual)
  }

  # As escalas variam entre painéis; a comparação principal é dentro de cada série.


  # Valores observados e previsões de um dia à frente

  curvas <- previsoes |>
    filter(ano == ano_grafico, localidade %in% localidades_grafico,
           familia %in% familias_grafico)

  observados <- curvas |>
    distinct(localidade, variavel, data_previsao, observado)

  grafico_previsao <- ggplot() +
    geom_line(data = observados,
              aes(data_previsao, observado, color = "Observado"), linewidth = 0.6) +
    geom_line(data = curvas,
              aes(data_previsao, previsto, color = familia), linewidth = 0.35) +
    facet_wrap(vars(localidade, variavel), ncol = 2, scales = "free_y") +
    scale_x_date(date_breaks = "3 months", date_labels = "%m/%Y",
                 expand = expansion(mult = 0)) +
    scale_color_manual(values = cores, breaks = c("Observado", familias_grafico)) +
    labs(x = NULL, y = "Umidade do solo", color = NULL,
         subtitle = paste("Previsões de um dia à frente em", ano_grafico)) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = linhas_legenda, byrow = TRUE))

  ggsave(file.path(pasta_figuras, paste0(comparacao, "_observado_previsto.pdf")),
         grafico_previsao, width = 10, height = 8)
  if (interactive()) print(grafico_previsao)

  # As curvas usam previsões já exportadas. Não há novos ajustes para este gráfico.
  # A sobreposição visual deve ser interpretada junto com RMSE e MAE.

}


# Resíduos de treinamento dos modelos selecionados

ajustes_ar <- readRDS("resultados/modelagem/autorregressivos/ajustes_selecionados.rds") |>
  mutate(familia = "AR")

ajustes_arma <- readRDS("resultados/modelagem/arma_arima/ajustes_selecionados.rds") |>
  mutate(familia = "ARMA/ARIMA")

ajustes_precipitacao <- readRDS(
  "resultados/modelagem/precipitacao/ajustes_selecionados.rds"
) |>
  semi_join(selecionados_precipitacao, by = chaves_modelo)

ajustes_setar <- readRDS(
  "resultados/modelagem/setar/ajustes_selecionados.rds"
) |>
  semi_join(selecionados_setar, by = chaves_modelo) |>
  mutate(familia = "SETAR")

ajustes_lstar <- readRDS(
  "resultados/modelagem/lstar/ajustes_selecionados.rds"
) |>
  semi_join(selecionados_lstar, by = chaves_modelo) |>
  mutate(familia = "LSTAR")

ajustes <- bind_rows(ajustes_ar, ajustes_arma, ajustes_precipitacao,
                     ajustes_setar, ajustes_lstar)
cadastro <- distinct(previsoes, across(all_of(chaves_serie)))
lista_residuos <- list()

# Os ajustes são lidos dos scripts 13 a 17. Nenhum modelo é reestimado aqui.

for (i in seq_len(nrow(ajustes))) {
  registro <- ajustes[i, ]

  lista_residuos[[i]] <- tibble(
    localidade_id = registro$localidade_id, variavel = registro$variavel,
    familia = registro$familia, modelo = registro$modelo,
    ano_validacao = registro$ano, fim_ajuste = registro$fim_ajuste,
    data = registro$dados_treino[[1]]$data,
    residuo = as.numeric(residuals(registro$ajuste[[1]]))
  )
}

residuos <- bind_rows(lista_residuos) |>
  left_join(cadastro, by = c("localidade_id", "variavel")) |>
  arrange(localidade_id, variavel, familia, ano_validacao, data)

# Cada ano de validação tem seu próprio ajuste e seus próprios resíduos.
# A persistência entra pela previsão; não possui um ajuste de parâmetros.

diagnosticos <- residuos |>
  group_by(across(all_of(chaves_serie)), familia, modelo, ano_validacao, fim_ajuste) |>
  group_modify(~ {
    valores_acf <- as.numeric(acf(.x$residuo, lag.max = 365, plot = FALSE)$acf)
    lag_pico <- which.max(abs(valores_acf[2:31]))
    tibble(n = nrow(.x), media = mean(.x$residuo), desvio_padrao = sd(.x$residuo),
           ACF_1 = valores_acf[2], ACF_7 = valores_acf[8],
           ACF_30 = valores_acf[31], ACF_365 = valores_acf[366],
           lag_pico_30 = lag_pico, ACF_pico_30 = valores_acf[lag_pico + 1L])
  }) |>
  ungroup()

print(diagnosticos, n = Inf, width = Inf)

# A ACF verifica a dependência que permaneceu após o ajuste.
# A média dos resíduos de treinamento não substitui o viés na validação.
# ACF_pico_30 é o maior coeficiente em módulo entre 1 e 30 dias, com seu sinal.
# No LSTAR, são os resíduos do modelo com a transição estimada no treinamento.

# No ajuste usado em 2024, o ARIMA de Belém/GWETTOP manteve ACF(2) = -0,132.
# Seu maior coeficiente residual em módulo foi maior que o do AR (0,081).
# Assim, o ganho de RMSE não significou menor dependência nos resíduos.

# Em Manaus/GWETTOP, o maior coeficiente em módulo caiu de 0,081 no AR
# para 0,038 no ARIMA. A melhora dos resíduos não foi uniforme entre séries.


# Nos ajustes usados em 2024, a interação apresentou maior pico de ACF
# residual em módulo, entre 1 e 30 dias, que o AR nas 12 séries.
# Em Petrolina/GWETTOP, esse pico passou de 0,025 para 0,110.
# Portanto, a redução do erro de previsão não eliminou a dependência residual.
# Os gráficos mostram resíduos de treinamento, não erros da validação.



# Nos ajustes usados em 2024, o SETAR apresentou menor pico de ACF residual
# em módulo, entre 1 e 30 dias, que o AR em oito das 12 séries.
# Em Belém/GWETTOP, esse pico caiu de 0,081 para 0,051.
# Em Petrolina/GWETTOP, aumentou de 0,025 para 0,038.


# Nos ajustes usados em 2024, o LSTAR apresentou menor pico de ACF residual
# em módulo, entre 1 e 30 dias, em sete séries frente ao AR e sete frente ao SETAR.
# Em Brasília/GWETTOP, o pico caiu de 0,039 no SETAR para 0,033 no LSTAR.
# Em Manaus/GWETTOP, aumentou de 0,077 para 0,081, apesar do menor RMSE.


# Gráficos de resíduos e ACF: uma página por série

residuos_grafico <- residuos |> filter(ano_validacao == ano_grafico)
series <- residuos_grafico |> group_by(localidade_id, variavel) |> group_split()

for (comparacao in names(comparacoes)) {
  familias_diagnostico <- setdiff(comparacoes[[comparacao]], "Persistência")

  pdf(file.path(pasta_diagnosticos, paste0("diagnosticos_", comparacao, "_residuos.pdf")),
      width = 10, height = max(8, 3 * length(familias_diagnostico)),
      pointsize = 10, onefile = TRUE)

  for (serie in series) {
    serie <- serie |> filter(familia %in% familias_diagnostico)
    par(mfrow = c(length(familias_diagnostico), 2),
        mar = c(4.2, 4.6, 3.1, 1), oma = c(0, 0, 2, 0), mgp = c(2.6, 0.8, 0))
    limite <- max(abs(serie$residuo))
    por_familia <- split(serie, serie$familia, drop = TRUE)
    acfs <- lapply(por_familia, function(dados) {
      as.numeric(acf(dados$residuo, lag.max = 30, plot = FALSE)$acf)[-1]
    })
    limite_acf <- 1.1 * max(0.05, abs(unlist(acfs)))

    for (familia_atual in familias_diagnostico) {
      dados_grafico <- por_familia[[familia_atual]]

      plot(dados_grafico$data, dados_grafico$residuo, type = "l",
           col = cores[[familia_atual]], ylim = c(-limite, limite),
           xlab = "", ylab = "Resíduo de treinamento", main = familia_atual)
      abline(h = 0, col = "grey50", lty = 2)

      plot(1:30, acfs[[familia_atual]], type = "h", col = cores[[familia_atual]],
           main = paste("ACF -", familia_atual), xlab = "Defasagem (dias)",
           ylab = "Autocorrelação", ylim = c(-limite_acf, limite_acf))
      abline(h = 0, col = "grey50")
      abline(h = c(-1, 1) * 1.96 / sqrt(nrow(dados_grafico)),
             col = "grey50", lty = 2)
    }

    mtext(paste(serie$localidade[1], "-", serie$variavel[1],
                "| Treinamento até", serie$fim_ajuste[1]), outer = TRUE, cex = 1)
  }
  dev.off()
}


# As linhas tracejadas da ACF são referências pontuais sob ruído branco.
# A defasagem zero foi omitida; a escala é a mesma para os modelos da página.
# Os gráficos detalham o ajuste usado em 2024; a tabela cobre os três ajustes anuais.

# Os resíduos apresentam picos e períodos com amplitudes diferentes.
# Essa inspeção, sozinha, não confirma variância não constante nem sua causa.


# Exportação para conferência e redação do artigo

write_csv(metricas, file.path(pasta_resultados, "metricas_validacao.csv"))
write_csv(metricas_anuais, file.path(pasta_resultados, "metricas_validacao_anuais.csv"))
write_csv(diagnosticos, file.path(pasta_resultados, "diagnosticos_residuos.csv"))

# Os resultados ainda pertencem à validação. O teste de 2025 permanece reservado.
