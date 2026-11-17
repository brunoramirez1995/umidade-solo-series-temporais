# ACF e PACF das Variáveis Alvo

library(tidyverse)
library(readr)
library(tidyr)
library(stringr)
library(purrr)

# Ler os dados

dados <- read_csv(
  "Dados/nasa_power_diario.csv",
  show_col_types = FALSE
)

dados$data <- as.Date(dados$data)


# Manter apenas o período de treinamento

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
  arrange(localidade, variavel, data)


# Ordem utilizada na análise

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


# Resumo de algumas autocorrelações importantes

calcular_autocorrelacoes <- function(local, alvo) {
  
  dados_serie <- alvos_treinamento |>
    filter(
      localidade == local,
      variavel == alvo
    ) |>
    arrange(data)
  
  valores_acf <- acf(
    dados_serie$umidade,
    lag.max = 366,
    plot = FALSE
  )$acf |>
    as.numeric()
  
  tibble(
    localidade = local,
    clima = first(dados_serie$clima),
    variavel = alvo,
    acf_1 = valores_acf[2],
    acf_7 = valores_acf[8],
    acf_30 = valores_acf[31],
    acf_365 = valores_acf[366],
    acf_366 = valores_acf[367]
  )
}


combinacoes <- expand_grid(
  localidade = localidades,
  variavel = variaveis
)


resumo_autocorrelacoes <- map2_dfr(
  combinacoes$localidade,
  combinacoes$variavel,
  calcular_autocorrelacoes
)


resumo_autocorrelacoes


# ACF e PACF das séries diárias

dir.create("figuras", showWarnings = FALSE)
dir.create("tabelas", showWarnings = FALSE)


pdf(
  "figuras/acf_pacf_series_diarias.pdf",
  width = 14,
  height = 5
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
    
    par(
      mfrow = c(1, 3),
      mar = c(4, 4, 3, 1),
      oma = c(0, 0, 2, 0)
    )
    
    # Dependência de curto prazo
    
    acf(
      serie,
      lag.max = 60,
      main = "ACF até 60 dias",
      xlab = "Defasagem (dias)"
    )
    
    # Dependência direta após controlar as defasagens intermediárias
    
    pacf(
      serie,
      lag.max = 60,
      main = "PACF até 60 dias",
      xlab = "Defasagem (dias)"
    )
    
    # Verificação do ciclo próximo de um ano
    
    acf(
      serie,
      lag.max = 400,
      main = "ACF até 400 dias",
      xlab = "Defasagem (dias)"
    )
    
    mtext(
      paste(local, "-", alvo),
      outer = TRUE,
      line = 0.5,
      font = 2,
      cex = 1.1
    )
  }
}


dev.off()

par(mfrow = c(1, 1))


# Salvar o resumo numérico

write_csv(
  resumo_autocorrelacoes,
  "tabelas/resumo_autocorrelacoes.csv"
)



# Todas as séries apresentam forte dependência temporal de curto prazo

# A autocorrelação de primeira ordem variou entre 0.958 e 0.999

# No atraso de 30 dias, a dependência foi maior nas localidades úmidas
# e menor nas localidades mais secas

# A ACF volta a apresentar valores positivos próximos dos atrasos 365 e 366,
# indicando repetição do padrão anual

# A PACF é dominada pelo primeiro atraso em todas as séries

# Em GWETROOT, o segundo atraso apresentou efeito negativo relevante,
# sugerindo que um AR(1) pode não representar toda a dependência de curto prazo

# Em GWETTOP, os atrasos seguintes variaram entre localidades,
# indicando possíveis diferenças na dinâmica temporal

# A ACF não apresenta corte abrupto, portanto um modelo MA puro
# não aparece como o principal candidato

# A queda lenta e oscilatória da ACF também é influenciada pela sazonalidade



# ACF e PACF selecionadas para o artigo
# Brasília, nas duas variáveis-alvo

pdf(
  "figuras/acf_pacf_brasilia_alvos.pdf",
  width = 10,
  height = 7.5
)

par(
  mfrow = c(2, 2),
  mar = c(4, 4, 3, 1),
  oma = c(0, 0, 2, 0)
)

for (alvo in variaveis) {
  
  serie_brasilia <- alvos_treinamento |>
    filter(
      localidade == "Brasília",
      variavel == alvo
    ) |>
    arrange(data) |>
    pull(umidade)
  
  acf(
    serie_brasilia,
    lag.max = 60,
    main = paste("ACF -", alvo),
    xlab = "Defasagem (dias)",
    ylim = c(0, 1)
  )
  
  pacf(
    serie_brasilia,
    lag.max = 60,
    main = paste("PACF -", alvo),
    xlab = "Defasagem (dias)",
    ylim = c(-0.5, 1)
  )
}

mtext(
  "Brasília - séries diárias",
  outer = TRUE,
  line = 0.5,
  font = 2,
  cex = 1.1
)

dev.off()

par(mfrow = c(1, 1))