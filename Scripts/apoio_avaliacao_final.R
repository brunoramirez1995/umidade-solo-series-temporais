# Funções de estimação utilizadas no script 20
# Mantêm as regras dos scripts 14, 16 e 17, sem executar suas grades.

# Reajusta uma configuração já escolhida, sem repetir a seleção.
# treino e avaliacao têm observado, sen_k, cos_k e tempo.

ajustar_prever_arma <- function(treino, avaliacao, config) {
  entradas <- c(
    paste0("sen_", seq_len(config$K)),
    paste0("cos_", seq_len(config$K)),
    if (config$tendencia) "tempo"
  )
  x_treino <- as.matrix(treino[, entradas, drop = FALSE])
  x_avaliacao <- as.matrix(avaliacao[, entradas, drop = FALSE])

  ajuste <- stats::arima(
    treino$observado, order = c(config$p, config$d, config$q),
    seasonal = c(0, 0, 0), xreg = x_treino,
    include.mean = config$d == 0, method = "ML",
    SSinit = "Rossignol2011", optim.control = list(maxit = 1000)
  )
  if (ajuste$code != 0L) stop("Ajuste sem convergência: código ", ajuste$code, ".")
  if (any(!is.finite(coef(ajuste)))) stop("Coeficientes não finitos.")
  if (any(!is.finite(ajuste$var.coef)) ||
      inherits(try(chol(ajuste$var.coef), silent = TRUE), "try-error")) {
    stop("Covariância dos coeficientes não positiva definida.")
  }

  parte_regressao <- as.numeric(
    x_avaliacao %*% coef(ajuste)[colnames(x_avaliacao)]
  )
  if ("intercept" %in% names(coef(ajuste))) {
    parte_regressao <- parte_regressao + coef(ajuste)["intercept"]
  }

  # Primeiro se prevê; só depois a observação atualiza o estado.
  # Os coeficientes ficam fixos durante todo o período avaliado.
  estado <- ajuste$model
  previsto <- numeric(nrow(avaliacao))
  for (t in seq_len(nrow(avaliacao))) {
    previsto[t] <- parte_regressao[t] +
      as.numeric(stats::KalmanForecast(1L, estado)$pred)
    atualizacao <- stats::KalmanRun(
      avaliacao$observado[t] - parte_regressao[t], estado,
      nit = -1L, update = TRUE
    )
    estado <- attr(atualizacao, "mod")
  }
  if (any(!is.finite(previsto))) stop("Previsões não finitas.")

  # Em d=1, arima() diferencia a relação completa, inclusive xreg.
  list(
    ajuste = ajuste, previsto = previsto,
    residuos = as.numeric(residuals(ajuste)), entradas = entradas
  )
}


# SETAR: cada regime tem intercepto e lags próprios; os demais termos são comuns.

matriz_setar <- function(base, config, limiar) {
  baixo <- as.numeric(base[[paste0("lag_", config$delta)]] <= limiar)
  alto <- 1 - baixo
  x_baixo <- cbind(1, as.matrix(base[paste0("lag_", seq_len(config$p1))])) * baixo
  x_alto <- cbind(1, as.matrix(base[paste0("lag_", seq_len(config$p2))])) * alto
  colnames(x_baixo) <- c("intercepto_baixo", paste0("phi_baixo_", seq_len(config$p1)))
  colnames(x_alto) <- c("intercepto_alto", paste0("phi_alto_", seq_len(config$p2)))
  comuns <- c(paste0("sen_", seq_len(config$K)), paste0("cos_", seq_len(config$K)),
              if (config$tendencia) "tempo")
  cbind(x_baixo, x_alto, as.matrix(base[comuns]))
}

ajustar_setar <- function(treino, config, fracao_minima = 0.15) {
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
  if (!is.finite(melhor_limiar)) stop("Nenhum limiar admissível com posto completo.")

  x <- matriz_setar(treino, config, melhor_limiar)
  ajuste <- lm(observado ~ . - 1, data = data.frame(observado = treino$observado, x))
  # A fórmula usa só as colunas do modelo; as matrizes temporárias não precisam ser salvas.
  attr(ajuste$terms, ".Environment") <- baseenv()
  attr(attr(ajuste$model, "terms"), ".Environment") <- baseenv()
  baixo <- z <= melhor_limiar
  list(ajuste = ajuste, limiar = melhor_limiar, soma_quadrados = melhor_soma,
       n_treino = length(z), n_baixo = sum(baixo), n_alto = sum(!baixo),
       transicoes_treino = sum(baixo[-1] != head(baixo, -1)))
}


# LSTAR: a logística combina as duas estruturas autorregressivas.

componentes_lstar <- function(base, config) {
  baixo <- cbind(1, as.matrix(base[paste0("lag_", seq_len(config$p1))]))
  alto <- cbind(1, as.matrix(base[paste0("lag_", seq_len(config$p2))]))
  colnames(baixo) <- c("intercepto_baixo", paste0("phi_baixo_", seq_len(config$p1)))
  colnames(alto) <- c("intercepto_alto", paste0("phi_alto_", seq_len(config$p2)))
  comuns <- c(paste0("sen_", seq_len(config$K)), paste0("cos_", seq_len(config$K)),
              if (config$tendencia) "tempo")
  list(baixo = baixo, alto = alto, comuns = as.matrix(base[comuns]))
}

matriz_lstar <- function(componentes, peso) {
  cbind(componentes$baixo * (1 - peso), componentes$alto * peso, componentes$comuns)
}

ajustar_lstar <- function(treino, config) {
  limites_gamma <- c(0.1, 100)
  gamas_iniciais <- c(1, 10, 50)
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

  # A regressão é resolvida em cada centro/inclinação; só os dois últimos são otimizados.
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

  centros <- unique(as.numeric(quantile(z_pad, c(0.25, 0.50, 0.75))))
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
                               message = conditionMessage(e)))
    avaliar(otim$par)
    tentativas[[j]] <- tibble::tibble(
      partida = j, centro_inicial = media_z + centros[j] * desvio_z,
      gamma_inicial = exp(inicial[1]), centro = media_z + otim$par[2] * desvio_z,
      gamma = exp(otim$par[1]), kappa = gamma / desvio_z,
      soma_quadrados = ultima_avaliacao$soma, convergencia = otim$convergence,
      posto_completo = ultima_avaliacao$valido,
      mensagem = if (is.null(otim$message)) "" else otim$message)
  }

  tentativas <- dplyr::bind_rows(tentativas)
  validas <- tentativas |>
    dplyr::filter(convergencia == 0, posto_completo, is.finite(soma_quadrados)) |>
    dplyr::arrange(soma_quadrados)
  if (nrow(validas) == 0) return(list(ajuste = NULL, tentativas = tentativas))
  melhor <- validas[1, ]
  peso <- plogis(melhor$kappa * (z - melhor$centro))
  x <- matriz_lstar(partes, peso)
  ajuste <- lm(observado ~ . - 1, data = data.frame(observado = y, x))
  attr(ajuste$terms, ".Environment") <- baseenv()
  attr(attr(ajuste$model, "terms"), ".Environment") <- baseenv()
  if (ajuste$rank != ncol(x)) stop("Coeficientes não identificados no ajuste final.")

  centro_pad <- (melhor$centro - media_z) / desvio_z
  list(ajuste = ajuste, tentativas = tentativas, centro = melhor$centro,
       gamma = melhor$gamma, kappa = melhor$kappa, media_transicao = media_z,
       desvio_transicao = desvio_z, soma_quadrados = melhor$soma_quadrados,
       n_convergentes = nrow(validas),
       limite_gamma = any(abs(log(melhor$gamma / limites_gamma)) < 1e-5),
       limite_centro = any(abs(centro_pad - limites_centro) < 1e-5),
       peso_min_treino = min(peso), peso_max_treino = max(peso),
       percentual_transicao_treino = 100 * mean(peso > 0.1 & peso < 0.9))
}

# Os parâmetros estimados ficam fixos; somente as entradas observadas são atualizadas.

prever_setar <- function(ajustado, base, config) {
  x <- matriz_setar(base, config, ajustado$limiar)
  as.numeric(predict(ajustado$ajuste, newdata = as.data.frame(x)))
}

prever_lstar <- function(ajustado, base, config) {
  z <- base[[paste0("lag_", config$delta)]]
  peso <- plogis(ajustado$kappa * (z - ajustado$centro))
  x <- matriz_lstar(componentes_lstar(base, config), peso)
  as.numeric(predict(ajustado$ajuste, newdata = as.data.frame(x)))
}
