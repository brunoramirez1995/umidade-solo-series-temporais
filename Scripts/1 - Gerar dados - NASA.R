
# Script - Gerar Dados

# Dados diários da NASA POWER

library(nasapower)
library(tidyverse)

# Pontos

(localidades <- tribble(
  ~localidade_id,     ~localidade,     ~uf, ~clima,               ~latitude, ~longitude,
  "manaus_am",        "Manaus",        "AM", "Úmido",                -3.1190,   -60.0217,
  "belem_pa",         "Belém",         "PA", "Úmido",                -1.4558,   -48.4902,
  "brasilia_df",      "Brasília",      "DF", "Tropical sazonal",    -15.7939,   -47.8828,
  "uberlandia_mg",    "Uberlândia",    "MG", "Tropical sazonal",    -18.9128,   -48.2755,
  "petrolina_pe",     "Petrolina",     "PE", "Mais seco",            -9.3891,   -40.5031,
  "montes_claros_mg", "Montes Claros", "MG", "Mais seco",           -16.7286,   -43.8582
))


# Alvo e covariáveis

parametros <- c(
  "GWETTOP",
  "GWETROOT",
  "PRECTOTCORR",
  "T2M",
  "RH2M",
  "ALLSKY_SFC_SW_DWN",
  "EVPTRNS"
)

# Download de um ponto

baixar_ponto <- function(localidade_id, localidade, uf, clima, latitude, longitude) {

  dados <- get_power(
    community = "ag",
    lonlat = c(longitude, latitude),
    pars = parametros,
    dates = c("2001-01-01", "2025-12-31"),
    temporal_api = "daily",
    time_standard = "LST"
  )

  Sys.sleep(1)

  dados |>
    as_tibble() |>
    transmute(
      localidade_id = localidade_id,
      localidade = localidade,
      uf = uf,
      clima = clima,
      latitude = latitude,
      longitude = longitude,
      data = YYYYMMDD,
      ano = YEAR,
      mes = MM,
      dia = DD,
      dia_ano = DOY,
      GWETTOP,
      GWETROOT,
      PRECTOTCORR,
      T2M,
      RH2M,
      ALLSKY_SFC_SW_DWN,
      EVPTRNS
    )
}

# Base final

dados_diarios <- pmap_dfr(localidades, baixar_ponto) |>
  arrange(localidade_id, data)

glimpse(dados_diarios)

# Dicionário das variáveis

(dicionario <- tribble(
  ~variavel,           ~unidade,    ~descricao,
  "GWETTOP",           "0 a 1",     "Saturação relativa do solo entre 0 e 5 cm",
  "GWETROOT",          "0 a 1",     "Saturação relativa do solo entre 0 e 100 cm",
  "PRECTOTCORR",       "mm/dia",    "Precipitação corrigida",
  "T2M",               "°C",         "Temperatura média a 2 m",
  "RH2M",              "%",          "Umidade relativa a 2 m",
  "ALLSKY_SFC_SW_DWN", "MJ/m²/dia", "Radiação solar incidente na superfície",
  "EVPTRNS",           "MJ/m²/dia", "Fluxo de energia de evapotranspiração"
))


# Salvar arquivos

dir.create("Dados", recursive = TRUE, showWarnings = FALSE)

write_csv(localidades, "dados/localidades.csv")

write_csv(dados_diarios, "dados/nasa_power_diario.csv")

write_csv(dicionario, "dados/dicionario_variaveis_diarias.csv")


