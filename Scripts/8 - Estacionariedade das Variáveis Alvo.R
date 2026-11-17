

# Estacionariedade das Variáveis Alvo

library(tidyverse)
library(readr)
library(tidyr)
library(stringr)
library(purrr)
library(tseries)



# Ler os dados

dados <- read_csv(
  "Dados/nasa_power_diario.csv",
  show_col_types = FALSE
)

dados$data <- as.Date(dados$data)


# Período de treinamento

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


# Ordem das séries

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



# Função para avaliar cada série

avaliar_estacionariedade <- function(local, alvo) {
  
  dados_serie <- alvos_treinamento |>
    filter(
      localidade == local,
      variavel == alvo
    ) |>
    arrange(data)
  
  clima_serie <- first(dados_serie$clima)
  serie_original <- dados_serie$umidade
  
  
  # Versões avaliadas
  
  series_avaliadas <- list(
    "Original" = serie_original,
    
    "Primeira diferença" = diff(
      serie_original,
      differences = 1
    ),
    
    "Diferença sazonal" = diff(
      serie_original,
      lag = 365
    ),
    
    "Primeira e sazonal" = diff(
      diff(
        serie_original,
        lag = 365
      ),
      differences = 1
    )
  )
  
  
  imap_dfr(
    series_avaliadas,
    function(serie, tipo) {
      
      # ADF: hipótese nula de raiz unitária
      
      teste_adf <- suppressWarnings(
        adf.test(
          serie,
          alternative = "stationary"
        )
      )
      
      
      # KPSS: hipótese nula de estacionariedade
      
      teste_kpss_nivel <- suppressWarnings(
        kpss.test(
          serie,
          null = "Level",
          lshort = FALSE
        )
      )
      
      teste_kpss_tendencia <- suppressWarnings(
        kpss.test(
          serie,
          null = "Trend",
          lshort = FALSE
        )
      )
      
      
      # Autocorrelação após cada diferenciação
      
      acf_1 <- acf(
        serie,
        lag.max = 1,
        plot = FALSE
      )$acf[2] |>
        as.numeric()
      
      
      tibble(
        localidade = local,
        clima = clima_serie,
        variavel = alvo,
        transformacao = tipo,
        observacoes = length(serie),
        media = mean(serie),
        desvio_padrao = sd(serie),
        acf_1 = acf_1,
        adf_estatistica = unname(teste_adf$statistic),
        adf_p_valor = teste_adf$p.value,
        kpss_nivel_estatistica =
          unname(teste_kpss_nivel$statistic),
        kpss_nivel_p_valor =
          teste_kpss_nivel$p.value,
        kpss_tendencia_estatistica =
          unname(teste_kpss_tendencia$statistic),
        kpss_tendencia_p_valor =
          teste_kpss_tendencia$p.value
      )
    }
  )
}




combinacoes <- expand_grid(
  localidade = localidades,
  variavel = variaveis
)


resultado_estacionariedade <- map2_dfr(
  combinacoes$localidade,
  combinacoes$variavel,
  avaliar_estacionariedade
)


# Organizar a apresentação

resultado_estacionariedade <- resultado_estacionariedade |>
  mutate(
    transformacao = factor(
      transformacao,
      levels = c(
        "Original",
        "Primeira diferença",
        "Diferença sazonal",
        "Primeira e sazonal"
      )
    )
  ) |>
  arrange(
    localidade,
    variavel,
    transformacao
  )


resultado_estacionariedade


resumo_estacionariedade <- resultado_estacionariedade |>
  mutate(
    resultado_adf = if_else(
      adf_p_valor < 0.05,
      "Rejeita raiz unitária",
      "Não rejeita raiz unitária"
    ),
    resultado_kpss_nivel = if_else(
      kpss_nivel_p_valor < 0.05,
      "Rejeita estacionariedade",
      "Não rejeita estacionariedade"
    ),
    resultado_kpss_tendencia = if_else(
      kpss_tendencia_p_valor < 0.05,
      "Rejeita estacionariedade",
      "Não rejeita estacionariedade"
    )
  ) |>
  select(
    localidade,
    clima,
    variavel,
    transformacao,
    acf_1,
    adf_p_valor,
    resultado_adf,
    kpss_nivel_p_valor,
    resultado_kpss_nivel,
    kpss_tendencia_p_valor,
    resultado_kpss_tendencia
  )


resumo_estacionariedade



dir.create(
  "tabelas",
  showWarnings = FALSE
)


write_csv(
  resultado_estacionariedade,
  "tabelas/resultados_testes_estacionariedade.csv"
)

write_csv(
  resumo_estacionariedade,
  "tabelas/resumo_estacionariedade.csv"
)




# O teste ADF rejeitou a hipótese de raiz unitária em todas as séries
# originais, indicando que a elevada persistência não caracteriza,
# isoladamente, um processo integrado.

# O KPSS de nível rejeitou a estacionariedade em 8 das 12 séries.
# Quando foi permitida uma tendência determinística, a hipótese de
# estacionariedade não foi rejeitada em 10 das 12 séries.

# Esses resultados sugerem que as mudanças na média estão mais
# relacionadas à tendência e à sazonalidade do que a uma raiz unitária.

# A primeira diferença produziu séries estacionárias segundo os testes,
# mas isso não significa que d = 1 deva ser adotado automaticamente.

# A diferença sazonal de período 365 manteve elevada autocorrelação no
# primeiro atraso e ainda apresentou rejeições no teste KPSS.

# A combinação das diferenças regular e sazonal não mostrou vantagem
# clara sobre a primeira diferença e reduz ainda mais o número de observações.

# Como decisão preliminar, serão comparados modelos com d = 0, após o
# tratamento da tendência e da sazonalidade, e modelos com d = 1.
# A diferenciação sazonal D = 1 não será adotada inicialmente.