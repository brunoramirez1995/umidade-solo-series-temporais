
# Verificação da Base

library(tidyverse)
library(readr)
library(tidyr)
library(stringr)
library(purrr)


# Ler arquivo

dados <- read_csv("Dados/nasa_power_diario.csv", show_col_types = FALSE)

glimpse(dados)

# Ajuste da data

dados$data <- as.Date(dados$data)

# Verificação da base

variaveis_numericas <- c(
  "GWETTOP",
  "GWETROOT",
  "PRECTOTCORR",
  "T2M",
  "RH2M",
  "ALLSKY_SFC_SW_DWN",
  "EVPTRNS"
)


# Chaves duplicadas

(duplicacoes <- dados |>
  count(localidade_id, data) |>
  filter(n > 1)
)

# Datas esperadas

(datas_esperadas <- expand_grid(
  localidade_id = unique(dados$localidade_id),
  data = seq(
    min(dados$data),
    max(dados$data),
    by = "day"
  )
)
)


(datas_ausentes <- datas_esperadas |>
  anti_join(
    dados |> distinct(localidade_id, data),
    by = c("localidade_id", "data")
  )
)


# Quantidade de observações por localidade

(contagem_localidades <- dados |>
  count(localidade_id)
)

# Cadastro das localidades

(
cadastro_localidades <- dados |>
  distinct(
    localidade_id,
    localidade,
    uf,
    clima,
    latitude,
    longitude
  )

)

# Resumo geral

verificacao_base <- tibble(
  observacoes = nrow(dados),
  localidades = n_distinct(dados$localidade_id),
  menor_serie = min(contagem_localidades$n),
  maior_serie = max(contagem_localidades$n),
  data_inicial = min(dados$data),
  data_final = max(dados$data),
  chaves_duplicadas = nrow(duplicacoes),
  datas_ausentes = nrow(datas_ausentes),
  valores_ausentes = sum(is.na(dados)),
  codigos_menos_999 = sum(
    as.matrix(
      dados |> select(all_of(variaveis_numericas))
    ) == -999,
    na.rm = TRUE
  ),
  valores_infinitos = sum(
    is.infinite(
      as.matrix(
        dados |> select(all_of(variaveis_numericas))
      )
    )
  ),
  inconsistencias_calendario = sum(
    dados$ano != as.integer(format(dados$data, "%Y")) |
      dados$mes != as.integer(format(dados$data, "%m")) |
      dados$dia != as.integer(format(dados$data, "%d")) |
      dados$dia_ano != as.integer(format(dados$data, "%j"))
  ),
  inconsistencias_cadastro =
    nrow(cadastro_localidades) -
    n_distinct(dados$localidade_id)
)

verificacao_base


# Verificação das variáveis numéricas

verificacao_variaveis <- dados |>
  select(all_of(variaveis_numericas)) |>
  pivot_longer(
    cols = everything(),
    names_to = "variavel",
    values_to = "valor"
  ) |>
  group_by(variavel) |>
  summarise(
    minimo = min(valor),
    maximo = max(valor),
    valores_distintos = n_distinct(valor),
    valores_ausentes = sum(is.na(valor)),
    codigos_menos_999 = sum(valor == -999, na.rm = TRUE),
    .groups = "drop"
  )

verificacao_variaveis


# Limites das variáveis alvo

limites_alvos <- dados |>
  select(localidade, GWETTOP, GWETROOT) |>
  pivot_longer(
    cols = c(GWETTOP, GWETROOT),
    names_to = "variavel",
    values_to = "umidade"
  ) |>
  group_by(localidade, variavel) |>
  summarise(
    minimo = min(umidade),
    maximo = max(umidade),
    iguais_a_zero = sum(umidade == 0),
    iguais_a_um = sum(umidade == 1),
    percentual_igual_um = mean(umidade == 1) * 100,
    fora_do_intervalo = sum(
      umidade < 0 | umidade > 1
    ),
    valores_distintos = n_distinct(umidade),
    .groups = "drop"
  )

limites_alvos



dir.create(
  "tabelas",
  showWarnings = FALSE,
  recursive = TRUE
)



write_csv(
  verificacao_base,
  "tabelas/verificacao_base.csv"
)

write_csv(
  verificacao_variaveis,
  "tabelas/verificacao_variaveis.csv"
)

write_csv(
  limites_alvos,
  "tabelas/limites_alvos.csv"
)