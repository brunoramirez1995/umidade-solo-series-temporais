
# Identificação Preliminar de Modelos ARMA e ARIMA

library(tidyverse)
library(readr)
library(tidyr)
library(stringr)
library(purrr)
library(tseries)



# Ler arquivo

dados <- read_csv(
  "Dados/nasa_power_diario.csv",
  show_col_types = FALSE
)

dados$data <- as.Date(dados$data)


# Utilizar somente o conjunto de treinamento

dados_treinamento <- dados |>
  filter(data <= as.Date("2021-12-31"))


# Organizar as variáveis alvo

alvos_treinamento <- dados_treinamento |>
  select(
    localidade,
    clima,
    data,
    GWETTOP,
    GWETROOT
  ) |>
  pivot_longer(
    cols = c(GWETTOP, GWETROOT),
    names_to = "variavel",
    values_to = "umidade"
  ) |>
  arrange(
    localidade,
    variavel,
    data
  )


# Ordem das localidades

localidades <- c(
  "Belém",
  "Manaus",
  "Brasília",
  "Uberlândia",
  "Montes Claros",
  "Petrolina"
)

variaveis <- c(
  "GWETTOP",
  "GWETROOT"
)

representacoes <- c(
  "Original",
  "Primeira diferença"
)


# Número de defasagens para a identificação de p e q

lag_maximo <- 30


# Calcular ACF e PACF

calcular_acf_pacf <- function(
    localidade_atual,
    variavel_atual,
    representacao_atual) {
  
  dados_serie <- alvos_treinamento |>
    filter(
      localidade == localidade_atual,
      variavel == variavel_atual
    ) |>
    arrange(data)
  
  serie <- dados_serie$umidade
  
  if (representacao_atual == "Primeira diferença") {
    serie <- diff(serie)
  }
  
  resultado_acf <- acf(
    serie,
    lag.max = lag_maximo,
    plot = FALSE
  )
  
  resultado_pacf <- pacf(
    serie,
    lag.max = lag_maximo,
    plot = FALSE
  )
  
  tibble(
    localidade = localidade_atual,
    clima = first(dados_serie$clima),
    variavel = variavel_atual,
    representacao = representacao_atual,
    observacoes = length(serie),
    defasagem = 1:lag_maximo,
    acf = as.numeric(resultado_acf$acf)[-1],
    pacf = as.numeric(resultado_pacf$acf),
    limite_aproximado = 1.96 / sqrt(length(serie))
  )
}


# Combinações analisadas

combinacoes <- expand_grid(
  localidade_atual = localidades,
  variavel_atual = variaveis,
  representacao_atual = representacoes
)


# Resultados para as 12 séries

(resultados_acf_pacf <- pmap_dfr(
  combinacoes,
  calcular_acf_pacf
))


# Resumo das primeiras nove defasagens

(resumo_acf_pacf <- resultados_acf_pacf |>
    filter(defasagem <= 9) |>
    select(
      localidade,
      clima,
      variavel,
      representacao,
      observacoes,
      defasagem,
      acf,
      pacf
    ) |>
    pivot_wider(
      names_from = defasagem,
      values_from = c(acf, pacf),
      names_sep = "_"
    ) |>
    arrange(
      variavel,
      localidade,
      representacao
    )
)


# Pastas

dir.create(
  "tabelas",
  showWarnings = FALSE
)

dir.create(
  "figuras",
  showWarnings = FALSE
)


# Exportar tabelas

write_csv(
  resultados_acf_pacf,
  "tabelas/acf_pacf_identificacao_completa.csv"
)

write_csv(
  resumo_acf_pacf,
  "tabelas/resumo_acf_pacf_baixa_ordem.csv"
)


# Gráficos completos para o apêndice

pdf(
  "figuras/acf_pacf_identificacao_baixa_ordem.pdf",
  width = 10,
  height = 8
)

for (local in localidades) {
  
  for (alvo in variaveis) {
    
    serie <- alvos_treinamento |>
      filter(
        localidade == local,
        variavel == alvo
      ) |>
      arrange(data) |>
      pull(umidade)
    
    serie_diferenciada <- diff(serie)
    
    par(
      mfrow = c(2, 2),
      mar = c(3, 3, 3, 1),
      oma = c(0, 0, 2, 0)
    )
    
    acf(
      serie,
      lag.max = lag_maximo,
      ylim = c(-0.6, 1),
      main = "ACF — série original"
    )
    
    pacf(
      serie,
      lag.max = lag_maximo,
      ylim = c(-0.6, 1),
      main = "PACF — série original"
    )
    
    acf(
      serie_diferenciada,
      lag.max = lag_maximo,
      ylim = c(-0.6, 1),
      main = "ACF — primeira diferença"
    )
    
    pacf(
      serie_diferenciada,
      lag.max = lag_maximo,
      ylim = c(-0.6, 1),
      main = "PACF — primeira diferença"
    )
    
    mtext(
      paste(local, "-", alvo),
      outer = TRUE,
      cex = 1.2,
      font = 2
    )
  }
}

dev.off()


# Figuras das localidades representativas

localidades_selecionadas <- c(
  "Belém",
  "Brasília",
  "Petrolina"
)

pdf(
  "figuras/acf_pacf_diferencas_selecionadas.pdf",
  width = 10,
  height = 9
)

for (alvo in variaveis) {
  
  par(
    mfrow = c(3, 2),
    mar = c(3, 3, 3, 1),
    oma = c(0, 0, 2, 0)
  )
  
  for (local in localidades_selecionadas) {
    
    serie_diferenciada <- alvos_treinamento |>
      filter(
        localidade == local,
        variavel == alvo
      ) |>
      arrange(data) |>
      pull(umidade) |>
      diff()
    
    acf(
      serie_diferenciada,
      lag.max = lag_maximo,
      ylim = c(-0.6, 0.6),
      xlim = c(0.5, lag_maximo + 0.5),
      xaxs = "i",
      xlab = "Defasagem (dias)",
      main = paste("ACF -", local)
    )
    
    pacf(
      serie_diferenciada,
      lag.max = lag_maximo,
      ylim = c(-0.6, 0.6),
      main = paste("PACF -", local)
    )
  }
  
  mtext(
    paste(alvo, "- primeira diferença"),
    outer = TRUE,
    cex = 1.2,
    font = 2
  )
}

dev.off()


# Nesta etapa não foi utilizado um procedimento automático
# para escolher as ordens p e q

# Em séries longas, correlações pequenas podem ultrapassar os
# limites de referência sem representar efeitos relevantes

# A identificação será baseada principalmente na forma da ACF
# e da PACF, considerando cortes e decaimentos nas primeiras defasagens

# Os candidatos serão mantidos com ordens baixas para preservar
# a interpretação e evitar modelos desnecessariamente complexos




# Nas séries originais, a ACF apresenta decaimento lento
# e a PACF é dominada pela primeira defasagem

# Esse comportamento confirma a forte persistência diária,
# mas ainda inclui os efeitos da tendência e da sazonalidade

# Após a primeira diferença, GWETTOP apresentou dependência
# mais curta e valores negativos frequentes na segunda defasagem

# GWETROOT manteve autocorrelações positivas mais elevadas,
# indicando maior memória temporal na camada radicular

# Não foi observado um mesmo corte da ACF ou da PACF
# em todas as localidades

# Por isso, foi definido um conjunto comum de modelos candidatos
# com ordens baixas, que será aplicado às 12 séries

# A seleção final não será baseada apenas nos gráficos
# e dependerá do ajuste, dos resíduos e da avaliação fora da amostra