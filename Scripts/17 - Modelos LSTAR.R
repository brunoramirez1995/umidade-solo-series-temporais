# 17 - Modelos LSTAR: transição suave entre dois regimes

# Pacotes e configurações

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
source("Scripts/apoio_modelagem.R", encoding = "UTF-8")

arquivo_dados <- "Dados/nasa_power_diario.csv"
pasta_resultados <- "resultados/modelagem/lstar"
pasta_figuras <- "figuras/principais"
limites_gamma <- c(0.1, 100)
gamas_iniciais <- c(1, 10, 50)
quantis_centros <- c(0.25, 0.50, 0.75)
inicio_grafico <- as.Date("2024-10-01")
fim_grafico <- as.Date("2024-12-31")
localidades_grafico <- c("Belém", "Brasília", "Petrolina")
chaves <- c("localidade_id", "localidade", "clima", "variavel", "modelo")
chaves_modelo <- c("localidade_id", "variavel", "modelo")
cores <- c("AR" = "#009E73", "SETAR" = "#7B3294",
           "LSTAR" = "#A65628", "Observado" = "grey25")

dir.create(pasta_resultados, recursive = TRUE, showWarnings = FALSE)
dir.create(pasta_figuras, recursive = TRUE, showWarnings = FALSE)


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
# Com delta = 1, a transição depende de Y_t. A chuva não entra neste modelo.
# A resposta permanece na escala original; a logística apenas define os pesos.


# Candidatos definidos na metodologia

configuracoes <- expand_grid(p1 = 1:3, p2 = 1:3, delta = 1:3,
                            K = 1:3, tendencia = c(FALSE, TRUE)) |>
  mutate(modelo = paste0("LSTAR", p1, "_", p2, "_D", delta, "_K", K,
                         if_else(tendencia, "_T", "")),
         n_parametros = 5L + p1 + p2 + 2L * K + as.integer(tendencia)) |>
  arrange(n_parametros, modelo)

# São 162 configurações por série, como no SETAR, com três ajustes anuais.
# A contagem inclui dois interceptos, centro, inclinação e variância dos erros.
# As ordens, os harmônicos e a tendência são selecionados novamente para o LSTAR.


# Componentes autorregressivos e termos comuns

componentes_lstar <- function(base, config) {
  baixo <- cbind(1, as.matrix(base[paste0("lag_", seq_len(config$p1))]))
  alto <- cbind(1, as.matrix(base[paste0("lag_", seq_len(config$p2))]))
  colnames(baixo) <- c("intercepto_baixo", paste0("phi_baixo_", seq_len(config$p1)))
  colnames(alto) <- c("intercepto_alto", paste0("phi_alto_", seq_len(config$p2)))
  comuns <- c(paste0("sen_", seq_len(config$K)),
              paste0("cos_", seq_len(config$K)), if (config$tendencia) "tempo")
  list(baixo = baixo, alto = alto, comuns = as.matrix(base[comuns]))
}

matriz_lstar <- function(componentes, peso) {
  cbind(componentes$baixo * (1 - peso), componentes$alto * peso,
        componentes$comuns)
}

# G próximo de zero favorece o componente baixo; próximo de um, o alto.
# Os termos harmônicos e a tendência têm coeficientes comuns aos componentes.


# Estimação da transição somente no treinamento

ajustar_lstar <- function(treino, config) {
  y <- treino$observado
  z <- treino[[paste0("lag_", config$delta)]]
  media_z <- mean(z)
  desvio_z <- sd(z)
  if (!is.finite(desvio_z) || desvio_z == 0) stop("Variável de transição constante.")
  z_pad <- (z - media_z) / desvio_z
  limites_centro <- as.numeric(quantile(z_pad, c(0.15, 0.85)))
  partes <- componentes_lstar(treino, config)
  n_baixo <- ncol(partes$baixo)
  n_alto <- ncol(partes$alto)
  ultimo_par <- NULL
  ultima_avaliacao <- NULL

  avaliar <- function(par, gradiente = FALSE) {
    if (!identical(par, ultimo_par)) {
      gamma <- exp(par[1])
      eta <- gamma * (z_pad - par[2])
      peso <- plogis(eta)
      x <- matriz_lstar(partes, peso)
      ajuste <- lm.fit(x, y)
      soma <- sum(ajuste$residuals^2)
      valido <- ajuste$rank == ncol(x) && is.finite(soma) && soma > 0 &&
        all(is.finite(ajuste$coefficients))
      valor <- 1e10
      derivadas <- c(0, 0)

      if (valido) {
        baixo <- drop(partes$baixo %*% ajuste$coefficients[seq_len(n_baixo)])
        alto <- drop(partes$alto %*% ajuste$coefficients[n_baixo + seq_len(n_alto)])
        auxiliar <- ajuste$residuals * (alto - baixo) * peso * (1 - peso)
        valor <- log(soma)
        derivadas <- c(-2 * sum(auxiliar * eta), 2 * gamma * sum(auxiliar)) / soma
      }

      ultimo_par <<- par
      ultima_avaliacao <<- list(valor = valor, gradiente = derivadas,
                                soma = soma, valido = valido)
    }
    if (gradiente) ultima_avaliacao$gradiente else ultima_avaliacao$valor
  }

  # O log da soma dos quadrados tem o mesmo mínimo e melhora a escala numérica.
  # O gradiente evita repetir regressões para aproximar as duas derivadas.
  # A última regressão é reutilizada quando optim pede valor e gradiente no mesmo ponto.

  centros <- unique(as.numeric(quantile(z_pad, quantis_centros)))
  tentativas <- vector("list", length(centros))

  for (j in seq_along(centros)) {
    valores <- vapply(gamas_iniciais, function(gamma) avaliar(c(log(gamma), centros[j])),
                      numeric(1))
    inicial <- c(log(gamas_iniciais[which.min(valores)]), centros[j])
    otim <- tryCatch(
      optim(inicial, avaliar, gr = function(par) avaliar(par, TRUE),
            method = "L-BFGS-B",
            lower = c(log(limites_gamma[1]), limites_centro[1]),
            upper = c(log(limites_gamma[2]), limites_centro[2]),
            control = list(maxit = 300, factr = 1e7, pgtol = 1e-8)),
      error = function(e) list(par = inicial, convergence = 99L,
                               message = conditionMessage(e))
    )
    avaliar(otim$par)
    tentativas[[j]] <- tibble(
      partida = j, centro_inicial = media_z + centros[j] * desvio_z,
      gamma_inicial = exp(inicial[1]),
      centro = media_z + otim$par[2] * desvio_z, gamma = exp(otim$par[1]),
      kappa = gamma / desvio_z, soma_quadrados = ultima_avaliacao$soma,
      convergencia = otim$convergence, posto_completo = ultima_avaliacao$valido,
      mensagem = if (is.null(otim$message)) "" else otim$message
    )
  }

  tentativas <- bind_rows(tentativas)
  validas <- tentativas |>
    filter(convergencia == 0, posto_completo, is.finite(soma_quadrados)) |>
    arrange(soma_quadrados)
  if (nrow(validas) == 0) return(list(ajuste = NULL, tentativas = tentativas))
  melhor <- validas[1, ]
  peso <- plogis(melhor$kappa * (z - melhor$centro))
  x <- matriz_lstar(partes, peso)
  ajuste <- lm(observado ~ . - 1, data = data.frame(observado = y, x))
  if (ajuste$rank != ncol(x)) stop("Coeficientes não identificados no ajuste final.")

  centro_pad <- (melhor$centro - media_z) / desvio_z
  list(ajuste = ajuste, tentativas = tentativas, centro = melhor$centro,
       gamma = melhor$gamma, kappa = melhor$kappa, desvio_transicao = desvio_z,
       n_convergentes = nrow(validas),
       limite_gamma = any(abs(log(melhor$gamma / limites_gamma)) < 1e-5),
       limite_centro = any(abs(centro_pad - limites_centro) < 1e-5),
       peso_min_treino = min(peso), peso_max_treino = max(peso),
       percentual_transicao_treino = 100 * mean(peso > 0.1 & peso < 0.9))
}

# São até três partidas, uma por centro inicial; repetições de quantis são removidas.
# Em cada centro, escolhe-se o melhor gamma inicial por erro de treinamento.
# gamma usa o desvio-padrão do treinamento; kappa = gamma / desvio preserva
# G = plogis(kappa * (z - centro)), como na equação do artigo.
# O centro fica entre os quantis 15% e 85%; isso não impõe ocupações mínimas dos pesos.
# Os limites de gamma são numéricos. Atingi-los é registrado, sem exclusão automática.


# Ajustes anuais e seleção na validação

series <- dados |> group_by(localidade_id, variavel) |> group_split()
lista_metricas <- list()
lista_anuais <- list()
lista_previsoes <- list()
lista_ajustes <- list()
lista_otimizacao <- list()
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
        ajustado <- ajustar_lstar(treino, config)
        lista_otimizacao[[length(lista_otimizacao) + 1L]] <- ajustado$tentativas |>
          mutate(localidade_id = serie$localidade_id[1], variavel = serie$variavel[1],
                 modelo = config$modelo, ano = janela$ano, .before = 1)
        if (is.null(ajustado$ajuste)) stop("Nenhuma partida convergiu com posto completo.")
        z <- avaliacao[[paste0("lag_", config$delta)]]
        peso <- plogis(ajustado$kappa * (z - ajustado$centro))
        x <- matriz_lstar(componentes_lstar(avaliacao, config), peso)
        previsto <- as.numeric(predict(ajustado$ajuste, newdata = as.data.frame(x)))
        if (any(!is.finite(previsto))) stop("Previsões não finitas.")

        list(
          previsoes = avaliacao |>
            transmute(localidade_id, localidade, clima, variavel,
                      data_origem = data - 1, data_previsao = data,
                      ano = janela$ano, familia = "LSTAR", modelo = config$modelo,
                      observado, previsto = .env$previsto, erro = observado - previsto,
                      transicao = .env$z, peso_alto = .env$peso,
                      centro = ajustado$centro, kappa = ajustado$kappa),
          ajuste = tibble(
            localidade_id = serie$localidade_id[1], localidade = serie$localidade[1],
            variavel = serie$variavel[1], familia = "LSTAR", modelo = config$modelo,
            p1 = config$p1, p2 = config$p2, delta = config$delta, K = config$K,
            tendencia = config$tendencia, ano = janela$ano,
            inicio_ajuste = janela$inicio_ajuste, fim_ajuste = janela$fim_ajuste,
            centro = ajustado$centro, gamma = ajustado$gamma, kappa = ajustado$kappa,
            desvio_transicao = ajustado$desvio_transicao,
            limite_gamma = ajustado$limite_gamma, limite_centro = ajustado$limite_centro,
            n_convergentes = ajustado$n_convergentes, n_treino = nrow(treino),
            peso_min_treino = ajustado$peso_min_treino,
            peso_max_treino = ajustado$peso_max_treino,
            percentual_transicao_treino = ajustado$percentual_transicao_treino,
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
    write_csv(falhas, file.path(pasta_resultados, "falhas.csv"))
    write_csv(bind_rows(lista_otimizacao), file.path(pasta_resultados, "otimizacao.csv"))
    stop(paste("Nenhuma configuração completou a validação para",
               serie$localidade[1], serie$variavel[1], "- confira falhas.csv."))
  }
  lista_previsoes[[length(lista_previsoes) + 1L]] <- melhor_previsao
  lista_ajustes[[length(lista_ajustes) + 1L]] <- melhores_ajustes
}

# O centro, a inclinação e os coeficientes são reestimados no início de cada ano.
# Dentro do ano, os parâmetros ficam fixos; os lags observados atualizam a previsão.
# A solução de cada ajuste é escolhida pelo treinamento; a configuração, pela validação.
# Apenas configurações com os três anos completos podem ser selecionadas.

metricas <- bind_rows(lista_metricas) |> left_join(configuracoes, by = "modelo")
metricas_anuais <- bind_rows(lista_anuais)
previsoes <- bind_rows(lista_previsoes)
ajustes_selecionados <- bind_rows(lista_ajustes)
otimizacao <- bind_rows(lista_otimizacao)
selecionados <- metricas |>
  filter(completo) |>
  group_by(localidade_id, variavel) |>
  arrange(RMSE, n_parametros, modelo, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup()


# Comparação com AR, SETAR e persistência

selecionados_ar <- read_csv(
  "resultados/modelagem/autorregressivos/selecionados_validacao.csv", show_col_types = FALSE)
selecionados_setar <- read_csv(
  "resultados/modelagem/setar/selecionados_validacao.csv", show_col_types = FALSE)
referencias <- bind_rows(
  selecionados_ar |> mutate(referencia = "AR"),
  selecionados_setar |> mutate(referencia = "SETAR"),
  read_csv("resultados/modelagem/referencias/metricas_validacao.csv", show_col_types = FALSE) |>
    filter(modelo == "persistencia") |> mutate(referencia = "persistencia")
) |>
  select(localidade_id, variavel, referencia, RMSE, MAE) |>
  pivot_wider(names_from = referencia, values_from = c(RMSE, MAE))

selecionados <- selecionados |>
  left_join(referencias, by = c("localidade_id", "variavel")) |>
  mutate(ganho_RMSE_AR = 100 * (1 - RMSE / RMSE_AR),
         ganho_MAE_AR = 100 * (1 - MAE / MAE_AR),
         ganho_RMSE_SETAR = 100 * (1 - RMSE / RMSE_SETAR),
         ganho_MAE_SETAR = 100 * (1 - MAE / MAE_SETAR),
         ganho_RMSE_persistencia = 100 * (1 - RMSE / RMSE_persistencia),
         ganho_MAE_persistencia = 100 * (1 - MAE / MAE_persistencia))

previsoes_ar <- read_csv(
  "resultados/modelagem/autorregressivos/previsoes_validacao.csv", show_col_types = FALSE) |>
  semi_join(selecionados_ar, by = chaves_modelo) |> mutate(familia = "AR")
previsoes_setar <- read_csv(
  "resultados/modelagem/setar/previsoes_validacao.csv", show_col_types = FALSE) |>
  semi_join(selecionados_setar, by = chaves_modelo) |> mutate(familia = "SETAR")

conferencia <- bind_rows(previsoes, previsoes_ar, previsoes_setar) |>
  group_by(localidade_id, variavel, data_previsao) |>
  summarise(n = n(), familias = n_distinct(familia), observados = n_distinct(observado),
            .groups = "drop")
if (any(conferencia$n != 3 | conferencia$familias != 3 | conferencia$observados != 1)) {
  stop("Datas ou observações não coincidem entre LSTAR, AR e SETAR.")
}


# Transição e coeficientes dos selecionados

pesos_validacao <- previsoes |>
  group_by(localidade_id, variavel, modelo, ano) |>
  summarise(peso_min_validacao = min(peso_alto), peso_max_validacao = max(peso_alto),
            percentual_transicao_validacao = 100 * mean(peso_alto > 0.1 & peso_alto < 0.9),
            .groups = "drop")
transicoes <- ajustes_selecionados |>
  select(-dados_treino, -ajuste) |>
  mutate(largura_10_90 = 2 * log(9) / kappa) |>
  left_join(pesos_validacao, by = c(chaves_modelo, "ano"))
coeficientes <- bind_rows(lapply(seq_len(nrow(ajustes_selecionados)), function(i) {
  registro <- ajustes_selecionados[i, ]
  valores <- coef(registro$ajuste[[1]])
  tibble(localidade_id = registro$localidade_id, localidade = registro$localidade,
         variavel = registro$variavel, modelo = registro$modelo, ano = registro$ano,
         termo = names(valores), coeficiente = unname(valores))
}))

print(selecionados |>
        select(localidade, variavel, p1, p2, delta, K, tendencia, RMSE, MAE,
               ganho_RMSE_AR, ganho_RMSE_SETAR, ganho_MAE_SETAR, vies, n_fora),
      n = Inf, width = Inf)
print(metricas_anuais |> semi_join(selecionados, by = chaves_modelo), n = Inf, width = Inf)
print(transicoes, n = Inf, width = Inf)
print(falhas, n = Inf, width = Inf)

# largura_10_90 mede quanto a umidade precisa variar para G passar de 0,1 a 0,9.
# Valores pequenos indicam transição mais abrupta; G é um peso, não uma probabilidade.
# Confira limites numéricos e tentativas com erros próximos e parâmetros diferentes.
# Os objetos lm são condicionados à transição estimada: seus testes usuais não
# representam a incerteza completa do LSTAR. Aqui são usados para previsão e resíduos.


# Gráficos de desempenho, previsão e transição

ganhos <- selecionados |>
  select(localidade, variavel, RMSE = ganho_RMSE_SETAR, MAE = ganho_MAE_SETAR) |>
  pivot_longer(c(RMSE, MAE), names_to = "metrica", values_to = "ganho") |>
  mutate(metrica = factor(metrica, levels = c("RMSE", "MAE")))
grafico_metricas <- ggplot(ganhos, aes(ganho, localidade)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5, color = cores[["LSTAR"]]) +
  facet_wrap(vars(metrica, variavel), ncol = 2, scales = "free_x") +
  labs(x = "Redução do erro em relação ao SETAR (%)", y = NULL) +
  theme_minimal(base_size = 11)

curvas <- bind_rows(previsoes, previsoes_ar, previsoes_setar) |>
  filter(data_previsao >= inicio_grafico, data_previsao <= fim_grafico,
         localidade %in% localidades_grafico)
observados <- curvas |> distinct(localidade, variavel, data_previsao, observado)
grafico_previsao <- ggplot() +
  geom_line(data = curvas, aes(data_previsao, previsto, color = familia), linewidth = 0.5) +
  geom_line(data = observados, aes(data_previsao, observado, color = "Observado"),
            linewidth = 0.4, linetype = "dashed") +
  facet_wrap(vars(localidade, variavel), ncol = 2, scales = "free_y") +
  scale_x_date(date_breaks = "1 month", date_labels = "%m/%Y", expand = expansion(mult = 0)) +
  scale_color_manual(values = cores, breaks = c("Observado", "AR", "SETAR", "LSTAR")) +
  labs(x = NULL, y = "Umidade do solo", color = NULL,
       subtitle = paste("Previsões de um dia à frente:", format(inicio_grafico, "%d/%m/%Y"),
                        "a", format(fim_grafico, "%d/%m/%Y"))) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

grafico_transicao <- previsoes |>
  filter(ano == 2024, localidade %in% localidades_grafico) |>
  ggplot(aes(transicao, peso_alto)) +
  geom_line(color = cores[["LSTAR"]], linewidth = 0.6) +
  geom_point(size = 0.6, color = cores[["LSTAR"]]) +
  geom_vline(data = transicoes |> filter(ano == 2024, localidade %in% localidades_grafico),
             aes(xintercept = centro), linetype = "dashed", color = "grey50") +
  facet_wrap(vars(localidade, variavel), ncol = 2, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = "Umidade usada na transição", y = "Peso do componente alto (G)",
       subtitle = "Valores da transição em 2024; parâmetros estimados até dezembro de 2023") +
  theme_minimal(base_size = 11)

grafico_pesos <- previsoes |>
  filter(ano == 2024, localidade %in% localidades_grafico) |>
  ggplot(aes(data_previsao, peso_alto)) +
  geom_line(color = cores[["LSTAR"]], linewidth = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  facet_wrap(vars(localidade, variavel), ncol = 2) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_x_date(date_breaks = "3 months", date_labels = "%m/%Y") +
  labs(x = "Dia previsto", y = "Peso do componente alto (G)",
       subtitle = "Pesos em 2024; parâmetros estimados até dezembro de 2023") +
  theme_minimal(base_size = 11)

# O recorte gráfico não altera as métricas dos três anos completos.
# A curva logística mostra o peso G, não a resposta prevista do solo.
if (interactive()) {
  print(grafico_metricas)
  print(grafico_previsao)
  print(grafico_transicao)
  print(grafico_pesos)
}


# O LSTAR reduziu o RMSE frente ao AR nas 12 séries e frente ao SETAR em oito.
# Os ganhos adicionais foram mais claros na camada superficial.
# Na zona radicular, as diferenças de RMSE em relação ao SETAR foram pequenas.
#
# A melhora no RMSE não foi acompanhada por redução do MAE em todos os casos.
# Algumas transições atingiram os limites da busca ou apresentaram pouca
# separação dos componentes. Seus parâmetros exigem cautela na interpretação.
#
# As previsões fora de [0,1] foram mantidas e registradas na avaliação.


# Exportação para os diagnósticos e o artigo

write_csv(metricas, file.path(pasta_resultados, "metricas_validacao.csv"))
write_csv(metricas_anuais, file.path(pasta_resultados, "metricas_validacao_anuais.csv"))
write_csv(selecionados, file.path(pasta_resultados, "selecionados_validacao.csv"))
write_csv(previsoes, file.path(pasta_resultados, "previsoes_validacao.csv"))
write_csv(transicoes, file.path(pasta_resultados, "transicoes_selecionadas.csv"))
write_csv(coeficientes, file.path(pasta_resultados, "coeficientes_selecionados.csv"))
write_csv(otimizacao, file.path(pasta_resultados, "otimizacao.csv"))
write_csv(falhas, file.path(pasta_resultados, "falhas.csv"))
saveRDS(ajustes_selecionados, file.path(pasta_resultados, "ajustes_selecionados.rds"))
ggsave(file.path(pasta_figuras, "lstar_validacao.pdf"), grafico_metricas, width = 9, height = 6.5)
ggsave(file.path(pasta_figuras, "lstar_observado_previsto.pdf"), grafico_previsao, width = 10, height = 8)
ggsave(file.path(pasta_figuras, "lstar_transicao_selecionadas.pdf"), grafico_transicao,
       width = 10, height = 8)
ggsave(file.path(pasta_figuras, "lstar_pesos_selecionadas.pdf"), grafico_pesos,
       width = 10, height = 8)

# Ganhos positivos indicam redução do erro frente à referência indicada na coluna.
# As previsões não são cortadas em [0,1]. O teste de 2025 permanece reservado.
