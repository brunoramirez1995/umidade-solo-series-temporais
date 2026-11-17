# 19 - Transferência entre localidades

# Pacotes e configurações

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

source("Scripts/apoio_modelagem.R", encoding = "UTF-8")

# Use "teste" para 2025 ou "validacao" para reproduzir 2022-2024.
# As duas etapas mantêm seus resultados em pastas separadas.

etapa <- "teste"
stopifnot(etapa %in% c("teste", "validacao"))

if (etapa == "teste") {
  calendario <- tibble(
    ano = 2025L, inicio_avaliacao = as.Date("2025-01-01"),
    fim_avaliacao = as.Date("2025-12-31")
  )
  ajustes <- readRDS("resultados/modelagem/teste_2025/ajustes_finais.rds") |>
    filter(familia == "AR + P + interação", ano == 2025L)
  pasta_resultados <- "resultados/modelagem/transferencia/teste_2025"
  periodo_completo <- "2025"
} else {
  calendario <- calendario_validacao
  ajustes <- readRDS("resultados/modelagem/precipitacao/ajustes_selecionados.rds") |>
    filter(com_interacao, ano %in% calendario$ano)
  pasta_resultados <- "resultados/modelagem/transferencia"
  periodo_completo <- "2022-2024"
}

pasta_figuras <- "figuras/principais"
pasta_previsoes <- "figuras/transferencia"
dir.create(pasta_resultados, recursive = TRUE, showWarnings = FALSE)
dir.create(pasta_figuras, recursive = TRUE, showWarnings = FALSE)
dir.create(pasta_previsoes, recursive = TRUE, showWarnings = FALSE)

inicio_grafico <- as.Date(paste0(max(calendario$ano), "-10-01"))
fim_grafico <- max(calendario$fim_avaliacao)
ordem_localidades <- c("Belém", "Manaus", "Brasília", "Uberlândia", "Montes Claros", "Petrolina")
tipos <- c("persistencia", "local", "direta", "calibrada")
rotulos <- c(persistencia = "Persistência", local = "Ajuste local",
             direta = "Transferência direta", calibrada = "Transferência calibrada")
cores <- c(persistencia = "#999999", local = "#009E73",
           direta = "#D55E00", calibrada = "#2166AC", Observado = "#222222")

# Os modelos locais e das origens já estão ajustados. Aqui são estimados
# somente o intercepto e a inclinação da calibração de cada transferência.
# Nenhuma origem ou configuração é escolhida com os resultados de 2025.


# Conferência dos ajustes salvos

stopifnot(
  nrow(ajustes) == 12L * nrow(calendario),
  !anyDuplicated(ajustes[c("localidade_id", "variavel", "ano")]),
  all(ajustes$inicio_ajuste == inicio_treinamento),
  all(ajustes$fim_ajuste == as.Date(paste0(ajustes$ano - 1, "-12-31")))
)
if (etapa == "teste") {
  stopifnot(all(vapply(ajustes$dados_treino, function(x) {
    max(x$data) <= as.Date("2024-12-31")
  }, logical(1))))
}


# Entradas do destino, construídas da mesma forma que no script 15

dados <- read_csv(
  "Dados/nasa_power_diario.csv",
  col_types = cols_only(
    localidade_id = col_character(), localidade = col_character(),
    clima = col_character(), data = col_date(),
    GWETTOP = col_double(), GWETROOT = col_double(), PRECTOTCORR = col_double()
  ),
  show_col_types = FALSE
) |>
  filter(data >= inicio_treinamento, data <= max(calendario$fim_avaliacao))

localidades <- dados |> distinct(localidade_id, localidade, clima)
stopifnot(nrow(localidades) == 6L, setequal(localidades$localidade_id, ajustes$localidade_id))

dados <- dados |>
  pivot_longer(all_of(alvos), names_to = "variavel", values_to = "observado") |>
  arrange(localidade_id, variavel, data) |>
  group_by(localidade_id, variavel) |>
  mutate(
    lag_1 = lag(observado, 1), lag_2 = lag(observado, 2), lag_3 = lag(observado, 3),
    precipitacao_lag_1 = lag(PRECTOTCORR, 1),
    interacao = lag_1 * precipitacao_lag_1
  ) |>
  ungroup() |>
  mutate(tempo = as.numeric(data - inicio_treinamento) / 365.25)

for (k in 1:3) {
  dados[[paste0("sen_", k)]] <- sin(2 * pi * k * dados$tempo)
  dados[[paste0("cos_", k)]] <- cos(2 * pi * k * dados$tempo)
}

# A linha é o dia previsto, t+1: lag_1 = Y_t e precipitacao_lag_1 = P_t.
# A data dos harmônicos e a referência de tempo são as mesmas da origem.
# Os valores observados atualizam as entradas, mas os coeficientes ficam fixos no ano.


# Aplicação direta e calibração no histórico do destino

lista_previsoes <- list()
lista_calibracoes <- list()
lista_dominio <- list()
falhas <- tibble(origem_id = character(), destino_id = character(),
                 variavel = character(), ano = integer(), motivo = character())

for (i in seq_len(nrow(ajustes))) {
  config <- ajustes[i, ]
  cadastro_origem <- localidades |> filter(localidade_id == config$localidade_id)
  destinos <- localidades |> filter(localidade_id != config$localidade_id)
  janela <- calendario |> filter(ano == config$ano)
  ajuste_origem <- config$ajuste[[1]]
  entradas <- all.vars(delete.response(terms(ajuste_origem)))

  message(cadastro_origem$localidade, " - ", config$variavel, " - ", config$ano)

  for (j in seq_len(nrow(destinos))) {
    cadastro_destino <- destinos[j, ]
    resultado <- tryCatch({
      ajuste_destino <- ajustes |>
        filter(localidade_id == cadastro_destino$localidade_id,
               variavel == config$variavel, ano == config$ano)
      if (nrow(ajuste_destino) != 1L) stop("Ajuste local ausente ou duplicado.")

      serie <- dados |>
        filter(localidade_id == cadastro_destino$localidade_id, variavel == config$variavel)
      treino <- serie |>
        filter(data >= config$inicio_ajuste, data <= config$fim_ajuste) |>
        drop_na(observado, all_of(entradas))
      avaliacao <- serie |>
        filter(data >= janela$inicio_avaliacao, data <= janela$fim_avaliacao)

      datas_esperadas <- seq(janela$inicio_avaliacao, janela$fim_avaliacao, by = "day")
      if (nrow(avaliacao) != length(datas_esperadas) ||
          anyNA(avaliacao$data) || !all(avaliacao$data == datas_esperadas)) {
        stop("Calendário do destino incompleto ou duplicado.")
      }

      # O modelo da origem mantém seus coeficientes e recebe as entradas do destino.
      treino$direta <- as.numeric(predict(ajuste_origem, newdata = treino))
      calibracao <- lm(observado ~ direta, data = treino)
      if (any(!is.finite(coef(calibracao)))) stop("Calibração não identificada.")

      previsoes <- avaliacao |>
        transmute(data_origem = data - 1, data_previsao = data, observado,
                  persistencia = lag_1,
                  local = as.numeric(predict(ajuste_destino$ajuste[[1]], newdata = avaliacao)),
                  direta = as.numeric(predict(ajuste_origem, newdata = avaliacao)))
      previsoes$calibrada <- as.numeric(predict(calibracao, newdata = previsoes))
      if (any(!is.finite(as.matrix(select(previsoes, observado, all_of(tipos)))))) {
        stop("Observações ou previsões não finitas.")
      }

      identificacao <- tibble(
        origem_id = cadastro_origem$localidade_id, origem = cadastro_origem$localidade,
        clima_origem = cadastro_origem$clima, destino_id = cadastro_destino$localidade_id,
        destino = cadastro_destino$localidade, clima_destino = cadastro_destino$clima,
        grupo = if_else(cadastro_origem$clima == cadastro_destino$clima,
                        "Mesmo grupo", "Outros grupos"),
        variavel = config$variavel, ano = config$ano
      )

      # As faixas vêm do treino da origem. São verificadas Y_t e P_t separadamente.
      dominio <- lapply(c("lag_1", "precipitacao_lag_1"), function(entrada) {
        x <- model.frame(ajuste_origem)[[entrada]]
        z <- treino[[entrada]]
        novo <- avaliacao[[entrada]]
        tibble(
          entrada = entrada, n_origem = length(x), n_destino = length(z),
          minimo_origem = min(x), q1_origem = unname(quantile(x, 0.25)),
          mediana_origem = median(x), q3_origem = unname(quantile(x, 0.75)),
          IQR_origem = IQR(x), maximo_origem = max(x),
          minimo_destino = min(z), q1_destino = unname(quantile(z, 0.25)),
          mediana_destino = median(z), q3_destino = unname(quantile(z, 0.75)),
          IQR_destino = IQR(z), maximo_destino = max(z),
          percentual_fora_historico = 100 * mean(z < min(x) | z > max(x)),
          percentual_fora_avaliacao = 100 * mean(novo < min(x) | novo > max(x))
        )
      }) |> bind_rows()

      list(
        previsoes = bind_cols(identificacao, previsoes),
        calibracao = bind_cols(identificacao, tibble(
          modelo_origem = config$modelo, modelo_destino = ajuste_destino$modelo,
          p_origem = config$p, K_origem = config$K, tendencia_origem = config$tendencia,
          p_destino = ajuste_destino$p, K_destino = ajuste_destino$K,
          tendencia_destino = ajuste_destino$tendencia,
          inicio_calibracao = min(treino$data), fim_calibracao = max(treino$data),
          n = nrow(treino), a = unname(coef(calibracao)[1]), b = unname(coef(calibracao)[2])
        )),
        dominio = bind_cols(identificacao, dominio)
      )
    }, error = function(e) e)

    if (inherits(resultado, "error")) {
      falhas <- bind_rows(falhas, tibble(
        origem_id = cadastro_origem$localidade_id, destino_id = cadastro_destino$localidade_id,
        variavel = config$variavel, ano = config$ano, motivo = conditionMessage(resultado)
      ))
    } else {
      lista_previsoes[[length(lista_previsoes) + 1L]] <- resultado$previsoes
      lista_calibracoes[[length(lista_calibracoes) + 1L]] <- resultado$calibracao
      lista_dominio[[length(lista_dominio) + 1L]] <- resultado$dominio
    }
  }
}

previsoes <- bind_rows(lista_previsoes)
calibracoes <- bind_rows(lista_calibracoes)
dominio <- bind_rows(lista_dominio)

write_csv(falhas, file.path(pasta_resultados, "falhas.csv"))
write_csv(previsoes, file.path(pasta_resultados, paste0("previsoes_", etapa, ".csv")))
write_csv(calibracoes, file.path(pasta_resultados, "calibracoes.csv"))
write_csv(dominio, file.path(pasta_resultados, "dominio.csv"))

# Uma falha não pode desaparecer na média nem retirar um par da comparação.
if (nrow(falhas) > 0L) stop("Confira falhas.csv antes de calcular os resumos.")
stopifnot(nrow(calibracoes) == 60L * nrow(calendario))

# A calibração usa respostas anteriores ao ano avaliado. Seus valores de ajuste
# não entram nas métricas. Não se impõe b positivo nem se limitam previsões a [0,1].


# Métricas por par, nas mesmas datas para as quatro alternativas

chaves_par <- c("origem_id", "origem", "clima_origem", "destino_id", "destino",
                "clima_destino", "grupo", "variavel")

resumir_previsoes <- function(grupos) {
  metricas <- previsoes |>
    pivot_longer(all_of(tipos), names_to = "tipo", values_to = "previsto") |>
    group_by(across(all_of(c(grupos, "tipo")))) |>
    summarise(calcular_metricas(observado, previsto),
              n_fora = sum(previsto < 0 | previsto > 1),
              percentual_fora = 100 * n_fora / n,
              excesso_maximo = max(pmax(-previsto, previsto - 1, 0)), .groups = "drop")

  referencias <- metricas |>
    filter(tipo %in% c("local", "persistencia")) |>
    select(all_of(grupos), tipo, RMSE, MAE) |>
    pivot_wider(names_from = tipo, values_from = c(RMSE, MAE))

  metricas |>
    left_join(referencias, by = grupos) |>
    mutate(ganho_RMSE_local = 100 * (1 - RMSE / RMSE_local),
           ganho_MAE_local = 100 * (1 - MAE / MAE_local),
           ganho_RMSE_persistencia = 100 * (1 - RMSE / RMSE_persistencia),
           ganho_MAE_persistencia = 100 * (1 - MAE / MAE_persistencia))
}

metricas <- resumir_previsoes(chaves_par)
metricas_anuais <- resumir_previsoes(c(chaves_par, "ano"))

# Ganho positivo: erro menor que o ajuste local da mesma família.
# O viés mantém a definição do projeto: observado - previsto.
# O local é repetido por par apenas para permitir a comparação nas mesmas datas.


# Comparação entre grupos climáticos, com o mesmo peso para cada destino

metricas_periodos <- metricas |> mutate(periodo = periodo_completo)
if (etapa == "validacao") {
  metricas_periodos <- bind_rows(
    metricas_periodos, metricas_anuais |> mutate(periodo = as.character(ano))
  )
}

ganhos_por_destino <- metricas_periodos |>
  filter(tipo %in% c("direta", "calibrada")) |>
  group_by(periodo, destino_id, destino, clima_destino, variavel, tipo, grupo) |>
  summarise(n_origens = n_distinct(origem_id),
            ganho_RMSE = mean(ganho_RMSE_local), ganho_MAE = mean(ganho_MAE_local),
            .groups = "drop")

stopifnot(all(ganhos_por_destino$n_origens ==
                ifelse(ganhos_por_destino$grupo == "Mesmo grupo", 1L, 4L)))

resumo_grupos <- ganhos_por_destino |>
  group_by(periodo, variavel, tipo, grupo) |>
  summarise(n_destinos = n(), ganho_RMSE_medio = mean(ganho_RMSE),
            ganho_RMSE_mediano = median(ganho_RMSE),
            ganho_MAE_medio = mean(ganho_MAE),
            destinos_com_ganho_RMSE = sum(ganho_RMSE > 0), .groups = "drop")

# Cada destino tem uma origem do mesmo grupo e quatro de outros grupos.
# Primeiro é calculada a média dessas quatro; depois, a média entre destinos.
# Os pares compartilham localidades e não são repetições independentes.

print(metricas |> filter(tipo %in% c("direta", "calibrada")) |>
        select(origem, destino, variavel, tipo, RMSE, MAE, ganho_RMSE_local,
               ganho_MAE_local, vies, n_fora), n = Inf, width = Inf)
print(resumo_grupos, n = Inf, width = Inf)




# Em 2025, as perdas médias foram menores entre localidades do mesmo grupo.
# Mesmo assim, todas as médias permaneceram negativas frente ao ajuste local.
# A semelhança climática não garantiu ganhos de transferência.


print(calibracoes |> select(origem, destino, variavel, ano, n, a, b),
      n = Inf, width = Inf)

# Na leitura dos resultados, confira se a calibração reduz o erro da aplicação
# direta e se chega a superar o modelo local. São duas comparações diferentes.
# Valores dentro das faixas de treino não garantem combinações semelhantes
# de entradas nem a mesma relação entre chuva e umidade.


# Resumos para interpretar os ganhos e as previsões fora dos limites

# Local e persistência aparecem cinco vezes nas métricas por par.
# Nesta contagem, cada uma das 12 séries locais entra uma única vez.
metricas_sem_repeticao <- bind_rows(
  metricas |> filter(tipo %in% c("direta", "calibrada")),
  metricas |> filter(tipo %in% c("local", "persistencia")) |>
    distinct(destino_id, variavel, tipo, .keep_all = TRUE)
)
resumo_alternativas <- metricas_sem_repeticao |>
  group_by(variavel, tipo) |>
  summarise(
    n_comparacoes = n(), n_previsoes = sum(n),
    ganhos_RMSE_local = sum(ganho_RMSE_local > 0),
    ganhos_MAE_local = sum(ganho_MAE_local > 0),
    ganhos_ambas_local = sum(ganho_RMSE_local > 0 & ganho_MAE_local > 0),
    ganhos_RMSE_persistencia = sum(ganho_RMSE_persistencia > 0),
    ganhos_MAE_persistencia = sum(ganho_MAE_persistencia > 0),
    ganhos_ambas_persistencia = sum(ganho_RMSE_persistencia > 0 & ganho_MAE_persistencia > 0),
    n_fora = sum(n_fora), percentual_fora = 100 * n_fora / n_previsoes,
    excesso_maximo = max(excesso_maximo), .groups = "drop"
  )

comparacao_calibracao <- metricas |>
  filter(tipo %in% c("direta", "calibrada")) |>
  select(all_of(chaves_par), tipo, RMSE, MAE, vies) |>
  pivot_wider(names_from = tipo, values_from = c(RMSE, MAE, vies)) |>
  mutate(ganho_RMSE = 100 * (1 - RMSE_calibrada / RMSE_direta),
         ganho_MAE = 100 * (1 - MAE_calibrada / MAE_direta),
         reducao_vies = abs(vies_calibrada) < abs(vies_direta))

resumo_calibracao <- comparacao_calibracao |>
  group_by(variavel) |>
  summarise(n_pares = n(), ganhos_RMSE = sum(ganho_RMSE > 0),
            ganhos_MAE = sum(ganho_MAE > 0), reducoes_vies = sum(reducao_vies),
            .groups = "drop")

print(resumo_alternativas, n = Inf, width = Inf)
print(resumo_calibracao, n = Inf, width = Inf)



# Resultados de 2025: a calibração reduziu o RMSE da aplicação direta
# em 52 das 60 combinações. Apenas 5 superaram o RMSE do modelo local.
# Melhorar a aplicação direta não significou superar o ajuste no destino.

# Em GWETROOT de Petrolina, três transferências calibradas superaram
# o modelo local no RMSE, mas todas permaneceram piores que a persistência.

# As contagens são descritivas: as origens e os destinos se repetem.





# Melhorar a transferência direta não implica superar o modelo local.
# As contagens dos pares são descritivas: as origens e os destinos se repetem.


# Gráficos dos ganhos por par e por grupo

ganhos <- metricas |>
  filter(tipo %in% c("direta", "calibrada")) |>
  mutate(origem = factor(origem, levels = rev(ordem_localidades)),
         destino = factor(destino, levels = ordem_localidades),
         tipo = factor(tipo, levels = c("direta", "calibrada"),
                       labels = rotulos[c("direta", "calibrada")]))

grafico_pares <- ggplot(ganhos, aes(destino, origem, fill = ganho_RMSE_local)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", ganho_RMSE_local)), size = 3) +
  facet_grid(tipo ~ variavel) +
  scale_fill_gradient2(low = "#E5A49B", mid = "white", high = "#8DB6D8", midpoint = 0) +
  labs(x = "Destino", y = "Origem", fill = "Ganho de RMSE (%)",
       subtitle = paste(periodo_completo, "- ganho em relação ao AR com precipitação e interação local")) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())

grafico_grupos <- ganhos_por_destino |>
  filter(periodo == periodo_completo) |>
  mutate(destino = factor(destino, levels = rev(ordem_localidades)),
         tipo = factor(tipo, levels = c("direta", "calibrada"))) |>
  ggplot(aes(ganho_RMSE, destino, color = grupo, shape = grupo)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  facet_wrap(vars(tipo, variavel), ncol = 2, scales = "free_x",
             labeller = labeller(tipo = rotulos)) +
  scale_color_manual(values = c("Mesmo grupo" = "#2166AC", "Outros grupos" = "#D55E00")) +
  labs(x = "Ganho de RMSE em relação ao ajuste local (%)", y = NULL, color = NULL,
       shape = NULL, subtitle = "Outros grupos: média das quatro origens em cada destino") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

ggsave(file.path(pasta_figuras, paste0("transferencia_ganhos_", etapa, ".pdf")),
       grafico_pares, width = 12, height = 9)
ggsave(file.path(pasta_figuras, paste0("transferencia_grupos_", etapa, ".pdf")),
       grafico_grupos, width = 11, height = 7)


# Observado e previsto: uma página por destino e camada, com as cinco origens

series_grafico <- previsoes |> distinct(destino_id, destino, variavel) |>
  arrange(match(destino, ordem_localidades), match(variavel, alvos))
pdf(file.path(pasta_previsoes, paste0("transferencia_observado_previsto_", etapa, ".pdf")),
    width = 11, height = 10, onefile = TRUE)

for (i in seq_len(nrow(series_grafico))) {
  serie <- series_grafico[i, ]
  curvas <- previsoes |>
    filter(destino_id == serie$destino_id, variavel == serie$variavel,
           data_previsao >= inicio_grafico, data_previsao <= fim_grafico)
  observados <- curvas |> select(origem, data_previsao, observado)
  curvas <- curvas |> pivot_longer(all_of(tipos), names_to = "tipo", values_to = "previsto")

  grafico <- ggplot(curvas, aes(data_previsao, previsto, color = tipo)) +
    geom_line(linewidth = 0.4) +
    geom_line(data = observados, aes(y = observado, color = "Observado"),
              linewidth = 0.45, linetype = "dashed") +
    facet_wrap(vars(origem), ncol = 1) +
    scale_color_manual(values = cores, breaks = c("Observado", tipos),
                       labels = c("Observado", unname(rotulos))) +
    scale_x_date(date_breaks = "1 month", date_labels = "%m/%Y",
                 expand = expansion(mult = 0)) +
    labs(title = paste("Destino:", serie$destino, "-", serie$variavel),
         subtitle = paste0("Outubro a dezembro de ", max(calendario$ano),
                           ". Cada painel corresponde a uma origem; horizonte de um dia."),
         x = NULL, y = "Umidade do solo", color = NULL) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  print(grafico)
}
dev.off()

# No teste, os gráficos mostram outubro a dezembro; as métricas usam todo 2025.
# Houve quatro previsões acima de um na aplicação direta e cinco na calibrada.
# São ocorrências por origem e destino, não necessariamente dias distintos.
# Esses valores foram mantidos na avaliação, sem truncamento.





# Exportação dos resumos para conferência e redação

write_csv(metricas, file.path(pasta_resultados, paste0("metricas_", etapa, ".csv")))
if (etapa == "validacao") {
  write_csv(metricas_anuais, file.path(pasta_resultados, "metricas_validacao_anuais.csv"))
}
write_csv(ganhos_por_destino, file.path(pasta_resultados, "ganhos_por_destino.csv"))
write_csv(resumo_grupos, file.path(pasta_resultados, "resumo_grupos_climaticos.csv"))
write_csv(resumo_alternativas, file.path(pasta_resultados, "resumo_alternativas.csv"))
write_csv(comparacao_calibracao, file.path(pasta_resultados, "comparacao_calibracao.csv"))
write_csv(resumo_calibracao, file.path(pasta_resultados, "resumo_calibracao.csv"))
capture.output(sessionInfo(), file = file.path(pasta_resultados, "sessionInfo.txt"))

# A comparação com a validação é acrescentada quando seu CSV já existe.
arquivo_validacao <- "resultados/modelagem/transferencia/metricas_validacao.csv"
if (etapa == "teste" && file.exists(arquivo_validacao)) {
  validacao <- read_csv(arquivo_validacao, show_col_types = FALSE) |>
    filter(tipo %in% c("direta", "calibrada")) |>
    select(origem_id, destino_id, variavel, tipo, RMSE, MAE,
           ganho_RMSE_local, ganho_MAE_local)
  comparacao_periodos <- metricas |>
    filter(tipo %in% c("direta", "calibrada")) |>
    select(all_of(chaves_par), tipo, RMSE, MAE, ganho_RMSE_local, ganho_MAE_local) |>
    left_join(validacao, by = c("origem_id", "destino_id", "variavel", "tipo"),
              suffix = c("_teste", "_validacao"))
  write_csv(comparacao_periodos, file.path(pasta_resultados, "comparacao_validacao_teste.csv"))
}

if (interactive()) {
  print(grafico_pares)
  print(grafico_grupos)
}

# Todas as direções são apresentadas, inclusive as que pioram o desempenho.
# O teste avalia o protocolo definido anteriormente; não altera suas escolhas.
