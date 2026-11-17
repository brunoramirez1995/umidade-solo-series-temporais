# 20 - Avaliação final dos modelos locais

# Pacotes e configurações

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

source("Scripts/apoio_avaliacao_final.R", encoding = "UTF-8")

inicio_treinamento <- as.Date("2001-01-01")
fim_treinamento <- as.Date("2024-12-31")
inicio_teste <- as.Date("2025-01-01")
fim_teste <- as.Date("2025-12-31")
datas_teste <- seq(inicio_teste, fim_teste, by = "day")
ano_teste <- as.integer(format(inicio_teste, "%Y"))

arquivo_configuracoes <- "resultados/modelagem/configuracoes_finais_2025.csv"
pasta_resultados <- "resultados/modelagem/teste_2025"
pasta_figuras <- "figuras/teste_2025"

dir.create(pasta_resultados, recursive = TRUE, showWarnings = FALSE)
dir.create(pasta_figuras, recursive = TRUE, showWarnings = FALSE)

inicio_grafico <- as.Date("2025-10-01")
fim_grafico <- fim_teste
alvos <- c("GWETTOP", "GWETROOT")
familias <- c("AR", "ARMA/ARIMA", "AR + P", "AR + P + interação", "SETAR", "LSTAR")
ordem_modelos <- c("Persistência", "Sazonal ingênua", familias)
ordem_localidades <- c("Belém", "Manaus", "Brasília", "Uberlândia", "Montes Claros", "Petrolina")
chaves <- c("localidade_id", "localidade", "clima", "variavel", "familia", "modelo")

# As configurações já foram escolhidas. Aqui são refeitos somente os 72 ajustes.
# Não há seleção pelo teste, mudança de família ou truncamento das previsões.


# Leitura das decisões e dos dados

configuracoes <- read_csv(arquivo_configuracoes, show_col_types = FALSE)
contagem <- configuracoes |> count(localidade_id, variavel)
if (nrow(configuracoes) != 72L || nrow(contagem) != 12L || any(contagem$n != 6L) ||
    !setequal(configuracoes$familia, familias) ||
    anyDuplicated(configuracoes[c("localidade_id", "variavel", "familia")])) {
  stop("Confira o arquivo congelado: são seis famílias para cada uma das 12 séries.")
}

dados <- read_csv(
  "Dados/nasa_power_diario.csv",
  col_types = cols_only(
    localidade_id = col_character(), localidade = col_character(),
    clima = col_character(), data = col_date(), GWETTOP = col_double(),
    GWETROOT = col_double(), PRECTOTCORR = col_double()
  ), show_col_types = FALSE
) |>
  filter(data >= inicio_treinamento, data <= fim_teste) |>
  arrange(localidade_id, data)

# Esta conferência protege os lags: uma linha deve corresponder a um dia.
calendario <- seq(inicio_treinamento, fim_teste, by = "day")
for (serie in split(dados, dados$localidade_id)) {
  if (!identical(serie$data, calendario) ||
      any(!is.finite(as.matrix(serie[c(alvos, "PRECTOTCORR")])))) {
    stop("Calendário incompleto ou valores ausentes em ", serie$localidade[1], ".")
  }
}
stopifnot(setequal(dados$localidade_id, configuracoes$localidade_id))

dados <- dados |>
  pivot_longer(all_of(alvos), names_to = "variavel", values_to = "observado") |>
  arrange(localidade_id, variavel, data) |>
  group_by(localidade_id, variavel) |>
  mutate(
    lag_1 = lag(observado, 1), lag_2 = lag(observado, 2), lag_3 = lag(observado, 3),
    sazonal_365 = lag(observado, 365), precipitacao_lag_1 = lag(PRECTOTCORR, 1),
    interacao = lag_1 * precipitacao_lag_1
  ) |>
  ungroup() |>
  mutate(tempo = as.numeric(data - inicio_treinamento) / 365.25)

for (k in 1:3) {
  dados[[paste0("sen_", k)]] <- sin(2 * pi * k * dados$tempo)
  dados[[paste0("cos_", k)]] <- cos(2 * pi * k * dados$tempo)
}

# Cada linha é o dia previsto: lag_1 = Y_t e precipitacao_lag_1 = P_t.
# Os lags são calculados antes do corte para preservar o histórico de dezembro.
# A sazonal repete 365 dias atrás; os harmônicos têm período de 365,25 dias.


# Persistência e referência sazonal

referencias <- dados |>
  filter(data %in% datas_teste) |>
  transmute(localidade_id, localidade, clima, variavel,
            data_origem = data - 1, data_previsao = data, ano = ano_teste,
            observado, persistencia = lag_1, sazonal_365) |>
  pivot_longer(c(persistencia, sazonal_365), names_to = "modelo", values_to = "previsto") |>
  mutate(familia = if_else(modelo == "persistencia", "Persistência", "Sazonal ingênua"),
         erro = observado - previsto)


# Reajuste até 2024 e previsões sucessivas de um dia

lista_previsoes <- list(referencias)
lista_ajustes <- list()
lista_residuos <- list()
lista_coeficientes <- list()
lista_otimizacao <- list()
falhas <- configuracoes[0, chaves] |> mutate(ano = integer(), motivo = character())

for (i in seq_len(nrow(configuracoes))) {
  config <- configuracoes[i, ]
  message(i, "/", nrow(configuracoes), " - ", config$localidade, " - ",
          config$variavel, " - ", config$familia)
  serie <- dados |>
    filter(localidade_id == config$localidade_id, variavel == config$variavel)
  treino <- serie |> filter(data <= fim_treinamento)
  avaliacao <- serie |> filter(data %in% datas_teste)

  resultado <- tryCatch({
    parametros <- tibble()
    transicao <- tibble()

    if (config$familia %in% c("AR", "AR + P", "AR + P + interação")) {
      entradas <- c(paste0("lag_", seq_len(config$p)),
                    paste0("sen_", seq_len(config$K)), paste0("cos_", seq_len(config$K)),
                    if (config$tendencia) "tempo",
                    if (config$com_precipitacao) "precipitacao_lag_1",
                    if (config$com_interacao) "interacao")
      treino <- treino |> drop_na(all_of(entradas))
      ajuste <- lm(reformulate(entradas, response = "observado"), data = treino)
      previsto <- as.numeric(predict(ajuste, newdata = avaliacao))

    } else if (config$familia == "ARMA/ARIMA") {
      # ARIMA recebe todo o histórico, sem descartar linhas por causa dos lags do AR.
      ajustado <- ajustar_prever_arma(treino, avaliacao, config)
      ajuste <- ajustado$ajuste
      previsto <- ajustado$previsto
      parametros <- tibble(convergencia = ajuste$code)

    } else {
      lags <- paste0("lag_", seq_len(max(config$p1, config$p2, config$delta)))
      treino <- treino |> drop_na(all_of(lags))
      z <- avaliacao[[paste0("lag_", config$delta)]]

      if (config$familia == "SETAR") {
        ajustado <- ajustar_setar(treino, config)
        previsto <- prever_setar(ajustado, avaliacao, config)
        transicao <- tibble(transicao = z, limiar = ajustado$limiar,
                            regime = if_else(z <= ajustado$limiar, "Baixo", "Alto"))
      } else {
        ajustado <- ajustar_lstar(treino, config)
        lista_otimizacao[[length(lista_otimizacao) + 1L]] <- bind_cols(
          config[rep(1L, nrow(ajustado$tentativas)), chaves], ajustado$tentativas)
        if (is.null(ajustado$ajuste)) stop("Nenhuma partida LSTAR convergiu com posto completo.")
        previsto <- prever_lstar(ajustado, avaliacao, config)
        transicao <- tibble(transicao = z, centro = ajustado$centro, kappa = ajustado$kappa,
                            peso_alto = plogis(ajustado$kappa * (z - ajustado$centro)))
      }
      ajuste <- ajustado$ajuste
      parametros <- as_tibble(ajustado[setdiff(names(ajustado), c("ajuste", "tentativas", "n_treino"))])
    }

    if (any(!is.finite(coef(ajuste)))) stop("Coeficientes não identificados.")
    if (length(previsto) != length(datas_teste) || any(!is.finite(previsto))) {
      stop("A configuração não produziu todas as previsões válidas.")
    }

    previsoes_modelo <- avaliacao |>
      transmute(localidade_id, localidade, clima, variavel,
                data_origem = data - 1, data_previsao = data, ano = ano_teste,
                familia = config$familia, modelo = config$modelo,
                observado, previsto = .env$previsto, erro = observado - previsto)
    if (ncol(transicao) > 0L) previsoes_modelo <- bind_cols(previsoes_modelo, transicao)

    registro <- config |>
      mutate(ano = ano_teste, inicio_ajuste = inicio_treinamento, fim_ajuste = fim_treinamento,
             n_treino = nrow(treino), dados_treino = list(select(treino, data, observado)),
             ajuste = list(ajuste))
    if (ncol(parametros) > 0L) registro <- bind_cols(registro, parametros)
    list(previsoes = previsoes_modelo, registro = registro,
         residuos = as.numeric(residuals(ajuste)), coeficientes = coef(ajuste))
  }, error = function(e) e)

  if (inherits(resultado, "error")) {
    falhas <- bind_rows(falhas, config[, chaves] |>
                         mutate(ano = ano_teste, motivo = conditionMessage(resultado)))
    message("  Falha registrada: ", conditionMessage(resultado))
    next
  }
  lista_previsoes[[length(lista_previsoes) + 1L]] <- resultado$previsoes
  lista_ajustes[[length(lista_ajustes) + 1L]] <- resultado$registro
  lista_residuos[[length(lista_residuos) + 1L]] <- bind_cols(
    config[rep(1L, nrow(treino)), chaves],
    tibble(data = treino$data, residuo = resultado$residuos))
  lista_coeficientes[[length(lista_coeficientes) + 1L]] <- bind_cols(
    config[rep(1L, length(resultado$coeficientes)), chaves],
    tibble(termo = names(resultado$coeficientes), estimativa = unname(resultado$coeficientes)))
}

# SETAR e LSTAR reestimam seus parâmetros de transição apenas no treinamento.
# Durante o teste, os lags observados atualizam as entradas, sem novo ajuste.
# ARIMA prevê primeiro e depois incorpora a observação ao filtro de Kalman.


# Salvar ajustes e resíduos para os scripts de transferência e diagnósticos

previsoes <- bind_rows(lista_previsoes)
ajustes_finais <- bind_rows(lista_ajustes)

write_csv(configuracoes, file.path(pasta_resultados, "configuracoes_utilizadas.csv"))
write_csv(previsoes, file.path(pasta_resultados, "previsoes_teste.csv"))
write_csv(falhas, file.path(pasta_resultados, "falhas.csv"))
saveRDS(bind_rows(lista_residuos), file.path(pasta_resultados, "residuos_treinamento.rds"))
write_csv(bind_rows(lista_coeficientes), file.path(pasta_resultados, "coeficientes.csv"))
write_csv(bind_rows(lista_otimizacao), file.path(pasta_resultados, "otimizacao_lstar.csv"))
saveRDS(ajustes_finais, file.path(pasta_resultados, "ajustes_finais.rds"))
if (nrow(ajustes_finais) > 0L) {
  write_csv(select(ajustes_finais, -dados_treino, -ajuste),
            file.path(pasta_resultados, "parametros_ajustes.csv"))
}
capture.output(sessionInfo(), file = file.path(pasta_resultados, "sessionInfo.txt"))

# Os objetos lm completos são preservados. O script 19 poderá reutilizar
# os 12 ajustes AR + P + interação, sem ajustá-los novamente.
# Os resíduos são do treinamento; erro em previsoes_teste.csv é do teste.


# Métricas nas mesmas datas, mantendo as falhas visíveis

esperados <- bind_rows(select(configuracoes, all_of(chaves)),
                       distinct(referencias, across(all_of(chaves))))
metricas <- previsoes |>
  group_by(across(all_of(chaves))) |>
  summarise(n = n(), RMSE = sqrt(mean(erro^2)), MAE = mean(abs(erro)), vies = mean(erro),
            n_fora = sum(previsto < 0 | previsto > 1),
            excesso_maximo = max(pmax(-previsto, previsto - 1, 0)),
            completo = identical(sort(data_previsao), datas_teste), .groups = "drop") |>
  right_join(esperados, by = chaves) |>
  mutate(n = coalesce(n, 0L), completo = coalesce(completo, FALSE),
         percentual_fora = if_else(completo, 100 * n_fora / n, NA_real_),
         across(c(RMSE, MAE, vies, excesso_maximo), ~ if_else(completo, .x, NA_real_)))

comparadores <- metricas |>
  filter(familia %in% c("Persistência", "AR")) |>
  mutate(referencia = if_else(familia == "AR", "AR", "persistencia")) |>
  select(localidade_id, variavel, referencia, RMSE, MAE) |>
  pivot_wider(names_from = referencia, values_from = c(RMSE, MAE))

# Se a referência tiver erro zero, o ganho percentual não é definido.
ganho_percentual <- function(erro, referencia) {
  ifelse(referencia > 0, 100 * (1 - erro / referencia), NA_real_)
}

metricas <- metricas |>
  left_join(comparadores, by = c("localidade_id", "variavel")) |>
  mutate(ganho_RMSE_persistencia = ganho_percentual(RMSE, RMSE_persistencia),
         ganho_MAE_persistencia = ganho_percentual(MAE, MAE_persistencia),
         ganho_RMSE_AR = ganho_percentual(RMSE, RMSE_AR),
         ganho_MAE_AR = ganho_percentual(MAE, MAE_AR)) |>
  arrange(match(localidade, ordem_localidades), match(variavel, alvos),
          match(familia, ordem_modelos))

write_csv(metricas, file.path(pasta_resultados, "metricas_teste.csv"))
print(select(metricas, localidade, variavel, familia, n, RMSE, MAE, vies,
             ganho_RMSE_persistencia, ganho_MAE_persistencia, n_fora, completo),
      n = Inf, width = Inf)
message(sum(metricas$completo), "/", nrow(metricas), " combinações completas; ",
        nrow(falhas), " falhas. Confira falhas.csv se houver registros.")

# Viés positivo indica subestimação média. Ganho positivo indica redução do erro.
# Uma configuração que falhe permanece na tabela, sem métricas na comparação principal.
# O teste descreve o desempenho das decisões anteriores; não escolhe novas configurações.


# Gráficos de comparação e previsão

# O mapa compara as seis famílias ajustadas com a persistência.
# As duas referências também constam da tabela e dos gráficos de previsão.

grafico_ganhos <- metricas |>
  filter(familia %in% familias) |>
  select(localidade, variavel, familia, completo,
         ganho_RMSE_persistencia, ganho_MAE_persistencia) |>
  pivot_longer(starts_with("ganho_"), names_to = "metrica", values_to = "ganho") |>
  mutate(metrica = factor(if_else(metrica == "ganho_RMSE_persistencia", "RMSE", "MAE"),
                         levels = c("RMSE", "MAE")),
         localidade = factor(localidade, levels = ordem_localidades),
         variavel = factor(variavel, levels = alvos),
         familia = factor(familia, levels = rev(familias)),
         rotulo = case_when(!completo ~ "Falha", is.finite(ganho) ~ sprintf("%.1f", ganho),
                             TRUE ~ "NA")) |>
  ggplot(aes(localidade, familia, fill = ganho)) +
  geom_tile(color = "white") + geom_text(aes(label = rotulo), size = 3) +
  facet_grid(metrica ~ variavel) +
  scale_fill_gradient2(low = "#D6604D", mid = "white", high = "#4393C3", midpoint = 0,
                       na.value = "grey85") +
  labs(x = NULL, y = NULL, fill = "Ganho (%)",
       title = paste(ano_teste, "- ganho em relação à persistência")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(pasta_figuras, "ganhos_teste.pdf"), grafico_ganhos, width = 12, height = 8)

# Uma página por série, com os oito modelos em painéis separados.
# O mesmo trimestre é mostrado em todas as séries, definido antes da avaliação.
pdf(file.path(pasta_figuras, "observado_previsto_teste.pdf"), width = 11, height = 10)
cadastro <- configuracoes |>
  distinct(localidade_id, localidade, variavel) |>
  arrange(match(localidade, ordem_localidades), match(variavel, alvos))
for (i in seq_len(nrow(cadastro))) {
  serie <- cadastro[i, ]
  quadro <- dados |>
    filter(localidade_id == serie$localidade_id, variavel == serie$variavel,
           data >= inicio_grafico, data <= fim_grafico) |>
    select(data_previsao = data, observado) |>
    crossing(familia = ordem_modelos) |>
    left_join(previsoes |>
                filter(localidade_id == serie$localidade_id, variavel == serie$variavel) |>
                select(data_previsao, familia, previsto), by = c("data_previsao", "familia")) |>
    mutate(familia = factor(familia, levels = ordem_modelos))
  falhas_serie <- metricas |>
    filter(localidade_id == serie$localidade_id, variavel == serie$variavel, !completo) |>
    mutate(familia = factor(familia, levels = ordem_modelos))
  grafico <- ggplot(quadro, aes(data_previsao)) +
    geom_line(aes(y = observado, color = "Observado"), linewidth = 0.45) +
    geom_line(aes(y = previsto, color = "Previsto"), linewidth = 0.4, na.rm = TRUE) +
    facet_wrap(~familia, ncol = 2) +
    scale_color_manual(values = c(Observado = "grey25", Previsto = "#2166AC")) +
    labs(x = NULL, y = "Umidade do solo", color = NULL,
         title = paste(serie$localidade, serie$variavel, ano_teste, sep = " - ")) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  if (nrow(falhas_serie) > 0L) {
    grafico <- grafico +
      geom_text(data = falhas_serie, aes(x = inicio_grafico, y = Inf),
                label = "Ajuste sem previsão completa", hjust = 0, vjust = 1.5,
                inherit.aes = FALSE, size = 3)
  }
  print(grafico)
}
dev.off()

# Os 72 ajustes produziram previsões completas para 2025, sem falhas.
# Na zona radicular, a precipitação reduziu os erros em cinco localidades.
# Petrolina foi a exceção: AR + P e AR + P + interação perderam para
# a persistência tanto no RMSE quanto no MAE.
#
# Na camada superficial, a interação reduziu o RMSE em relação ao
# AR + P nas seis localidades. Na zona radicular, isso ocorreu em três.
#
# Os ganhos de RMSE não foram acompanhados sempre por ganhos de MAE.
# A persistência permaneceu competitiva nas séries com muitos dias
# sem mudança no valor registrado.
#
# Em GWETROOT de Petrolina, dois episódios concentraram 47,5% da soma
# dos erros quadráticos do modelo com interação.
#
# Houve uma previsão acima de 1 e nenhuma previsão negativa.
# Os resultados de 2025 serão relatados sem modificar as configurações.
