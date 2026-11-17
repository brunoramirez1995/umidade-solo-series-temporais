# 16 - Modelos SETAR: dois regimes de umidade

# Pacotes e configurações

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)


source("Scripts/apoio_modelagem.R", encoding = "UTF-8")

arquivo_dados <- "Dados/nasa_power_diario.csv"
pasta_resultados <- "resultados/modelagem/setar"
pasta_figuras <- "figuras/principais"
fracao_minima <- 0.15
inicio_grafico <- as.Date("2024-10-01")
fim_grafico <- as.Date("2024-12-31")
localidades_grafico <- c("Belém", "Brasília", "Petrolina")
chaves <- c("localidade_id", "localidade", "clima", "variavel", "modelo")
chaves_modelo <- c("localidade_id", "variavel", "modelo")
cores <- c("AR" = "#009E73", "SETAR" = "#7B3294", "Observado" = "grey25")


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
  mutate(lag_1 = lag(observado, 1), lag_2 = lag(observado, 2),
         lag_3 = lag(observado, 3)) |>
  ungroup() |>
  mutate(tempo = as.numeric(data - inicio_treinamento) / 365.25)

for (k in 1:3) {
  dados[[paste0("sen_", k)]] <- sin(2 * pi * k * dados$tempo)
  dados[[paste0("cos_", k)]] <- cos(2 * pi * k * dados$tempo)
}

# Cada linha corresponde ao dia previsto, t+1.
# lag_delta é Y_(t+1-delta); com delta = 1, o regime depende de Y_t.
# Os harmônicos usam a data prevista. A precipitação não entra neste script.


# Candidatos definidos na metodologia

configuracoes <- expand_grid(p1 = 1:3, p2 = 1:3, delta = 1:3,
                            K = 1:3, tendencia = c(FALSE, TRUE)) |>
  mutate(modelo = paste0("SETAR", p1, "_", p2, "_D", delta, "_K", K,
                         if_else(tendencia, "_T", "")),
         n_parametros = 4L + p1 + p2 + 2L * K + as.integer(tendencia)) |>
  arrange(n_parametros, modelo)

# São 162 configurações por série e três ajustes anuais por configuração.
# A contagem inclui dois interceptos, o limiar e a variância dos erros.
# A ordenação resolve empates de RMSE pelo menor número de parâmetros.


# Construção do modelo: dinâmica própria e sazonalidade comum aos regimes

matriz_setar <- function(base, config, limiar) {
  baixo <- as.numeric(base[[paste0("lag_", config$delta)]] <= limiar)
  alto <- 1 - baixo

  x_baixo <- cbind(1, as.matrix(base[paste0("lag_", seq_len(config$p1))])) * baixo
  x_alto <- cbind(1, as.matrix(base[paste0("lag_", seq_len(config$p2))])) * alto
  colnames(x_baixo) <- c("intercepto_baixo", paste0("phi_baixo_", seq_len(config$p1)))
  colnames(x_alto) <- c("intercepto_alto", paste0("phi_alto_", seq_len(config$p2)))

  comuns <- c(paste0("sen_", seq_len(config$K)),
              paste0("cos_", seq_len(config$K)), if (config$tendencia) "tempo")
  cbind(x_baixo, x_alto, as.matrix(base[comuns]))
}

# Os indicadores ligam somente o intercepto e os lags do regime correspondente.
# Senos, cossenos e tendência entram uma única vez, com coeficientes comuns.


# Busca do limiar dentro do treinamento

ajustar_setar <- function(treino, config) {
  z <- treino[[paste0("lag_", config$delta)]]
  limiares <- sort(unique(z))
  n_baixo <- findInterval(limiares, sort(z))
  minimo <- ceiling(fracao_minima * length(z))
  limiares <- limiares[n_baixo >= minimo & length(z) - n_baixo >= minimo]

  melhor_soma <- Inf
  melhor_limiar <- NA_real_

  for (limiar in limiares) {
    x <- matriz_setar(treino, config, limiar)
    ajuste <- lm.fit(x, treino$observado)
    if (ajuste$rank < ncol(x)) next
    soma <- sum(ajuste$residuals^2)

    if (is.finite(soma) && soma < melhor_soma) {
      melhor_soma <- soma
      melhor_limiar <- limiar
    }
  }

  if (!is.finite(melhor_limiar)) stop("Nenhum limiar admissível com coeficientes identificados.")

  x <- matriz_setar(treino, config, melhor_limiar)
  base_ajuste <- data.frame(observado = treino$observado, x)
  ajuste <- lm(observado ~ . - 1, data = base_ajuste)
  list(ajuste = ajuste, limiar = melhor_limiar)
}

# A busca usa valores distintos: cortes entre dois valores consecutivos
# produzem a mesma divisão. Valores repetidos permanecem no mesmo regime.
# A restrição de 15% usa as contagens reais, e o critério é a soma dos
# quadrados dos resíduos no treinamento. A validação não escolhe o limiar.
# Não há intercepto global, pois cada regime já tem seu próprio intercepto.


# Ajustes anuais e seleção na validação

series <- dados |> group_by(localidade_id, variavel) |> group_split()
lista_metricas <- list()
lista_anuais <- list()
lista_previsoes <- list()
lista_ajustes <- list()
falhas <- tibble(localidade_id = character(), variavel = character(),
                 modelo = character(), ano = integer(), motivo = character())

for (serie in series) {
  message(serie$localidade[1], " - ", serie$variavel[1])
  melhor_RMSE <- Inf
  melhor_previsao <- NULL
  melhores_ajustes <- NULL

  for (i in seq_len(nrow(configuracoes))) {
    config <- configuracoes[i, ]
    lags <- paste0("lag_", seq_len(max(config$p1, config$p2, config$delta)))
    previsoes_config <- list()
    ajustes_config <- list()

    for (j in seq_len(nrow(calendario_validacao))) {
      janela <- calendario_validacao[j, ]
      treino <- serie |>
        filter(data >= janela$inicio_ajuste, data <= janela$fim_ajuste) |>
        drop_na(observado, all_of(lags))
      avaliacao <- serie |>
        filter(data >= janela$inicio_avaliacao, data <= janela$fim_avaliacao)

      resultado <- tryCatch({
        ajustado <- ajustar_setar(treino, config)
        x <- matriz_setar(avaliacao, config, ajustado$limiar)
        previsto <- as.numeric(predict(ajustado$ajuste, newdata = as.data.frame(x)))
        if (any(!is.finite(previsto))) stop("Previsões não finitas.")

        z <- treino[[paste0("lag_", config$delta)]]
        regime_treino <- ifelse(z <= ajustado$limiar, "Baixo", "Alto")

        list(
          previsoes = avaliacao |>
            transmute(localidade_id, localidade, clima, variavel,
                      data_origem = data - 1, data_previsao = data,
                      ano = janela$ano, familia = "SETAR", modelo = config$modelo,
                      observado, previsto = .env$previsto, erro = observado - previsto,
                      transicao = .data[[paste0("lag_", config$delta)]],
                      limiar = ajustado$limiar,
                      regime = if_else(transicao <= limiar, "Baixo", "Alto")),
          ajuste = tibble(
            localidade_id = serie$localidade_id[1], localidade = serie$localidade[1],
            variavel = serie$variavel[1], familia = "SETAR", modelo = config$modelo,
            p1 = config$p1, p2 = config$p2, delta = config$delta, K = config$K,
            tendencia = config$tendencia, ano = janela$ano,
            inicio_ajuste = janela$inicio_ajuste, fim_ajuste = janela$fim_ajuste,
            limiar = ajustado$limiar, n_treino = nrow(treino),
            n_baixo = sum(regime_treino == "Baixo"), n_alto = sum(regime_treino == "Alto"),
            transicoes = sum(regime_treino[-1] != head(regime_treino, -1)),
            dados_treino = list(select(treino, data, observado)),
            ajuste = list(ajustado$ajuste)
          )
        )
      }, error = function(e) e)

      if (inherits(resultado, "error")) {
        falhas <- bind_rows(falhas, tibble(
          localidade_id = serie$localidade_id[1], variavel = serie$variavel[1],
          modelo = config$modelo, ano = janela$ano, motivo = conditionMessage(resultado)))
      } else {
        previsoes_config[[j]] <- resultado$previsoes
        ajustes_config[[j]] <- resultado$ajuste
      }
    }

    if (length(previsoes_config) == 0) next
    previsoes_config <- bind_rows(previsoes_config)
    resumo <- previsoes_config |>
      group_by(across(all_of(chaves))) |>
      summarise(calcular_metricas(observado, previsto),
                n_fora = sum(previsto < 0 | previsto > 1),
                excesso_maximo = max(pmax(-previsto, previsto - 1, 0)), .groups = "drop") |>
      mutate(completo = n == as.integer(fim_validacao - inicio_validacao) + 1L,
             percentual_fora = 100 * n_fora / n)

    lista_metricas[[length(lista_metricas) + 1L]] <- resumo
    lista_anuais[[length(lista_anuais) + 1L]] <- previsoes_config |>
      group_by(across(all_of(c(chaves, "ano")))) |>
      summarise(calcular_metricas(observado, previsto),
                n_fora = sum(previsto < 0 | previsto > 1),
                excesso_maximo = max(pmax(-previsto, previsto - 1, 0)), .groups = "drop")

    if (resumo$completo && resumo$RMSE < melhor_RMSE) {
      melhor_RMSE <- resumo$RMSE
      melhor_previsao <- previsoes_config
      melhores_ajustes <- bind_rows(ajustes_config)
    }
    if (i %% 27 == 0) message("  ", i, " de ", nrow(configuracoes), " configurações")
  }

  if (is.null(melhor_previsao)) {
    dir.create(pasta_resultados, recursive = TRUE, showWarnings = FALSE)
    write_csv(falhas, file.path(pasta_resultados, "falhas.csv"))
    stop(paste("Nenhuma configuração completou a validação para",
               serie$localidade[1], serie$variavel[1], "- confira falhas.csv."))
  }

  lista_previsoes[[length(lista_previsoes) + 1L]] <- melhor_previsao
  lista_ajustes[[length(lista_ajustes) + 1L]] <- melhores_ajustes
}

# O limiar e os coeficientes são reestimados a cada ano, somente no treinamento.
# Dentro do ano, eles ficam fixos; os lags observados atualizam a previsão e o regime.
# As métricas de todos os candidatos são mantidas. Para economizar memória e espaço,
# apenas as previsões e os objetos dos selecionados são guardados, sem novos ajustes.

metricas <- bind_rows(lista_metricas) |> left_join(configuracoes, by = "modelo")
metricas_anuais <- bind_rows(lista_anuais)
previsoes <- bind_rows(lista_previsoes)
ajustes_selecionados <- bind_rows(lista_ajustes)

selecionados <- metricas |>
  filter(completo) |>
  group_by(localidade_id, variavel) |>
  arrange(RMSE, n_parametros, modelo, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup()


# Comparação com o AR e a persistência

configuracoes_ar <- read_csv(
  "resultados/modelagem/autorregressivos/selecionados_validacao.csv", show_col_types = FALSE)
referencias <- bind_rows(
  configuracoes_ar |> mutate(referencia = "AR"),
  read_csv("resultados/modelagem/referencias/metricas_validacao.csv",
           show_col_types = FALSE) |>
    filter(modelo == "persistencia") |> mutate(referencia = "persistencia")
) |>
  select(localidade_id, variavel, referencia, RMSE, MAE) |>
  pivot_wider(names_from = referencia, values_from = c(RMSE, MAE))

selecionados <- selecionados |>
  left_join(referencias, by = c("localidade_id", "variavel")) |>
  mutate(ganho_RMSE_AR = 100 * (1 - RMSE / RMSE_AR),
         ganho_MAE_AR = 100 * (1 - MAE / MAE_AR),
         ganho_RMSE_persistencia = 100 * (1 - RMSE / RMSE_persistencia),
         ganho_MAE_persistencia = 100 * (1 - MAE / MAE_persistencia))

previsoes_ar <- read_csv(
  "resultados/modelagem/autorregressivos/previsoes_validacao.csv", show_col_types = FALSE) |>
  semi_join(configuracoes_ar, by = chaves_modelo) |> mutate(familia = "AR")

conferencia <- previsoes |>
  left_join(previsoes_ar |>
              select(localidade_id, variavel, data_previsao, observado_AR = observado),
            by = c("localidade_id", "variavel", "data_previsao"))
if (anyNA(conferencia$observado_AR) || any(conferencia$observado != conferencia$observado_AR)) {
  stop("Datas ou observações não coincidem com a referência AR.")
}


# Limiares, ocupação dos regimes e coeficientes dos selecionados

regimes_validacao <- previsoes |>
  group_by(across(all_of(c(chaves, "ano")) )) |>
  summarise(n_validacao = n(), percentual_baixo_validacao = 100 * mean(regime == "Baixo"),
            transicoes_validacao = sum(regime[-1] != head(regime, -1)), .groups = "drop")

regimes <- ajustes_selecionados |>
  select(-dados_treino, -ajuste) |>
  mutate(percentual_baixo_treino = 100 * n_baixo / n_treino) |>
  left_join(regimes_validacao |> select(-localidade, -clima),
            by = c(chaves_modelo, "ano"))

coeficientes <- bind_rows(lapply(seq_len(nrow(ajustes_selecionados)), function(i) {
  registro <- ajustes_selecionados[i, ]
  valores <- coef(registro$ajuste[[1]])
  tibble(localidade_id = registro$localidade_id, localidade = registro$localidade,
         variavel = registro$variavel, modelo = registro$modelo, ano = registro$ano,
         termo = names(valores), coeficiente = unname(valores))
}))

print(selecionados |>
        select(localidade, variavel, p1, p2, delta, K, tendencia, RMSE, MAE,
               ganho_RMSE_AR, ganho_MAE_AR, ganho_RMSE_persistencia, vies, n_fora),
      n = Inf, width = Inf)

# O SETAR reduziu o RMSE do AR nas 12 séries e o MAE em oito.
# Os maiores ganhos ocorreram em GWETTOP de Brasília e Petrolina.
# Em GWETROOT, a redução do RMSE foi inferior a 0,8% nas seis localidades.

# Contra a persistência, houve melhora de RMSE nas 12 séries.
# Pelo MAE, a vantagem ocorreu apenas em Petrolina/GWETTOP.
# Dias sem mudança registrada favorecem a persistência nessa métrica.

# Foram obtidas 50 previsões acima de um, todas em GWETROOT:
# quatro em Belém e 46 em Manaus. O maior valor foi 1,00555.




print(metricas_anuais |> semi_join(selecionados, by = chaves_modelo), n = Inf, width = Inf)





# O RMSE foi menor que o do AR em 30 dos 36 casos anuais.
# Brasília e Petrolina/GWETTOP melhoraram nos três anos.
# Em Montes Claros/GWETTOP, houve pequena piora em 2024.



print(regimes, n = Inf, width = Inf)

# O limiar permaneceu igual nos três ajustes em oito das 12 séries.
# Os regimes acompanham parte do ciclo anual e não devem ser
# interpretados automaticamente como classes de solo seco e úmido.




print(falhas, n = Inf, width = Inf)

# Ganho positivo indica redução do erro. Confira RMSE, MAE e os três anos.
# Baixo e alto identificam os lados do limiar estimado, não classes físicas de solo.
# A restrição de 15% vale no treinamento; a ocupação pode mudar na validação.
# Limiares semelhantes entre anos não bastam para demonstrar estabilidade global.
# A seleção por previsão não constitui, sozinha, um teste de não linearidade.


# Gráficos de comparação, previsão e regimes

ganhos <- selecionados |>
  select(localidade, variavel, RMSE = ganho_RMSE_AR, MAE = ganho_MAE_AR) |>
  pivot_longer(c(RMSE, MAE), names_to = "metrica", values_to = "ganho") |>
  mutate(metrica = factor(metrica, levels = c("RMSE", "MAE")))
grafico_metricas <- ggplot(ganhos, aes(ganho, localidade)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5, color = cores[["SETAR"]]) +
  facet_wrap(vars(metrica, variavel), ncol = 2, scales = "free_x") +
  labs(x = "Redução do erro em relação ao AR (%)", y = NULL) + theme_minimal(base_size = 11)

curvas <- bind_rows(previsoes, previsoes_ar) |>
  filter(data_previsao >= inicio_grafico, data_previsao <= fim_grafico,
         localidade %in% localidades_grafico)
observados <- curvas |> distinct(localidade, variavel, data_previsao, observado)
grafico_previsao <- ggplot() +
  geom_line(data = curvas, aes(data_previsao, previsto, color = familia), linewidth = 0.5) +
  geom_line(data = observados, aes(data_previsao, observado, color = "Observado"),
            linewidth = 0.4, linetype = "dashed") +
  facet_wrap(vars(localidade, variavel), ncol = 2, scales = "free_y") +
  scale_x_date(date_breaks = "1 month", date_labels = "%m/%Y", expand = expansion(mult = 0)) +
  scale_color_manual(values = cores, breaks = c("Observado", "AR", "SETAR")) +
  labs(x = NULL, y = "Umidade do solo", color = NULL,
       subtitle = paste("Previsões de um dia à frente:", format(inicio_grafico, "%d/%m/%Y"),
                        "a", format(fim_grafico, "%d/%m/%Y"))) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

estados <- previsoes |>
  filter(ano == 2024, localidade %in% localidades_grafico)
limiares_grafico <- estados |> distinct(localidade, variavel, limiar)
grafico_regimes <- ggplot(estados, aes(data_previsao, transicao)) +
  geom_line(color = "grey65", linewidth = 0.3) +
  geom_point(aes(color = regime), size = 0.6) +
  geom_hline(data = limiares_grafico, aes(yintercept = limiar), linetype = "dashed") +
  facet_wrap(vars(localidade, variavel), ncol = 2, scales = "free_y") +
  scale_color_manual(values = c("Baixo" = "#D55E00", "Alto" = "#2166AC")) +
  scale_x_date(date_breaks = "3 months", date_labels = "%m/%Y", expand = expansion(mult = 0)) +
  labs(x = "Dia previsto", y = "Umidade usada para definir o regime", color = "Regime",
       subtitle = "Validação de 2024: limiares estimados até dezembro de 2023") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

# A figura dos regimes usa a variável de transição, não a umidade futura observada.
# O recorte das curvas não altera as métricas, calculadas nos três anos completos.

if (interactive()) {
  print(grafico_metricas)
  print(grafico_previsao)
  print(grafico_regimes)
}


# Exportação para os diagnósticos e o artigo

dir.create(pasta_resultados, recursive = TRUE, showWarnings = FALSE)
dir.create(pasta_figuras, recursive = TRUE, showWarnings = FALSE)
write_csv(metricas, file.path(pasta_resultados, "metricas_validacao.csv"))
write_csv(metricas_anuais, file.path(pasta_resultados, "metricas_validacao_anuais.csv"))
write_csv(selecionados, file.path(pasta_resultados, "selecionados_validacao.csv"))
write_csv(previsoes, file.path(pasta_resultados, "previsoes_validacao.csv"))
write_csv(regimes, file.path(pasta_resultados, "regimes_selecionados.csv"))
write_csv(coeficientes, file.path(pasta_resultados, "coeficientes_selecionados.csv"))
write_csv(falhas, file.path(pasta_resultados, "falhas.csv"))
saveRDS(ajustes_selecionados, file.path(pasta_resultados, "ajustes_selecionados.rds"))
ggsave(file.path(pasta_figuras, "setar_validacao.pdf"), grafico_metricas, width = 9, height = 6.5)
ggsave(file.path(pasta_figuras, "setar_observado_previsto.pdf"), grafico_previsao, width = 10, height = 8)
ggsave(file.path(pasta_figuras, "setar_regimes_selecionadas.pdf"), grafico_regimes, width = 10, height = 8)

# As configurações são selecionadas na validação. O teste de 2025 permanece reservado.
