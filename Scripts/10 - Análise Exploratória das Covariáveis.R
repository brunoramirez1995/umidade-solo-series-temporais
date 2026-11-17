# Análise Exploratória das Covariáveis

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

glimpse(dados)


# Conjunto de treinamento

dados_treinamento <- dados |>
  filter(data <= as.Date("2021-12-31"))


# Variáveis analisadas

covariaveis <- c(
  "PRECTOTCORR",
  "T2M",
  "RH2M",
  "ALLSKY_SFC_SW_DWN",
  "EVPTRNS"
)

alvos <- c(
  "GWETTOP",
  "GWETROOT"
)


# Ordem das localidades

ordem_localidades <- c(
  "Belém",
  "Manaus",
  "Brasília",
  "Uberlândia",
  "Montes Claros",
  "Petrolina"
)

ordem_climas <- c(
  "Úmido",
  "Tropical sazonal",
  "Mais seco"
)

localidades_selecionadas <- c(
  "Belém",
  "Brasília",
  "Petrolina"
)


# Nomes utilizados nos gráficos

rotulos_covariaveis <- c(
  PRECTOTCORR = "Precipitação (mm/dia)",
  T2M = "Temperatura a 2 m (°C)",
  RH2M = "Umidade relativa (%)",
  ALLSKY_SFC_SW_DWN = "Radiação solar (MJ/m²/dia)",
  EVPTRNS = "Evapotranspiração (MJ/m²/dia)"
)


# Organizar as covariáveis

covariaveis_treinamento <- dados_treinamento |>
  select(
    localidade_id,
    localidade,
    clima,
    data,
    ano,
    mes,
    all_of(covariaveis)
  ) |>
  pivot_longer(
    cols = all_of(covariaveis),
    names_to = "covariavel",
    values_to = "valor"
  ) |>
  mutate(
    localidade = factor(
      localidade,
      levels = ordem_localidades
    ),
    clima = factor(
      clima,
      levels = ordem_climas
    ),
    covariavel = factor(
      covariavel,
      levels = covariaveis
    )
  )


# Estatísticas descritivas

(resumo_covariaveis <- covariaveis_treinamento |>
    group_by(
      localidade,
      clima,
      covariavel
    ) |>
    summarise(
      observacoes = n(),
      media = mean(valor),
      desvio_padrao = sd(valor),
      minimo = min(valor),
      percentil_05 = unname(
        quantile(valor, 0.05)
      ),
      mediana = median(valor),
      percentil_95 = unname(
        quantile(valor, 0.95)
      ),
      maximo = max(valor),
      .groups = "drop"
    )
)


# Resumo específico da precipitação

(resumo_precipitacao <- dados_treinamento |>
    group_by(localidade, clima) |>
    summarise(
      observacoes = n(),
      precipitacao_media = mean(PRECTOTCORR),
      dias_sem_precipitacao = sum(
        PRECTOTCORR == 0
      ),
      percentual_sem_precipitacao = mean(
        PRECTOTCORR == 0
      ) * 100,
      dias_acima_1_mm = sum(
        PRECTOTCORR > 1
      ),
      percentual_acima_1_mm = mean(
        PRECTOTCORR > 1
      ) * 100,
      media_dias_com_precipitacao = mean(
        PRECTOTCORR[PRECTOTCORR > 0]
      ),
      maximo = max(PRECTOTCORR),
      .groups = "drop"
    )
)


# Relação das covariáveis com os alvos

dados_relacoes <- dados_treinamento |>
  arrange(localidade_id, data) |>
  group_by(localidade_id) |>
  mutate(
    variacao_atual_gwettop =
      GWETTOP - lag(GWETTOP),
    
    variacao_atual_gwetroot =
      GWETROOT - lag(GWETROOT),
    
    variacao_futura_gwettop =
      lead(GWETTOP) - GWETTOP,
    
    variacao_futura_gwetroot =
      lead(GWETROOT) - GWETROOT
  ) |>
  ungroup() |>
  pivot_longer(
    cols = all_of(covariaveis),
    names_to = "covariavel",
    values_to = "valor_covariavel"
  )


relacoes_gwettop <- dados_relacoes |>
  transmute(
    localidade,
    clima,
    data,
    alvo = "GWETTOP",
    covariavel,
    valor_covariavel,
    nivel_alvo = GWETTOP,
    variacao_atual = variacao_atual_gwettop,
    variacao_futura = variacao_futura_gwettop
  )


relacoes_gwetroot <- dados_relacoes |>
  transmute(
    localidade,
    clima,
    data,
    alvo = "GWETROOT",
    covariavel,
    valor_covariavel,
    nivel_alvo = GWETROOT,
    variacao_atual = variacao_atual_gwetroot,
    variacao_futura = variacao_futura_gwetroot
  )


relacoes_covariaveis_alvos <- bind_rows(
  relacoes_gwettop,
  relacoes_gwetroot
)


# Correlações por localidade

(correlacoes_covariaveis_alvos <-
    relacoes_covariaveis_alvos |>
    group_by(
      localidade,
      clima,
      alvo,
      covariavel
    ) |>
    summarise(
      correlacao_nivel = cor(
        valor_covariavel,
        nivel_alvo,
        use = "complete.obs"
      ),
      correlacao_variacao_atual = cor(
        valor_covariavel,
        variacao_atual,
        use = "complete.obs"
      ),
      correlacao_variacao_futura = cor(
        valor_covariavel,
        variacao_futura,
        use = "complete.obs"
      ),
      .groups = "drop"
    )
)


# Correlações entre as próprias covariáveis

pares_covariaveis <- combn(
  covariaveis,
  2,
  simplify = FALSE
)


(correlacoes_entre_covariaveis <-
    dados_treinamento |>
    group_by(localidade, clima) |>
    group_modify(
      ~ map_dfr(
        pares_covariaveis,
        function(par) {
          
          tibble(
            covariavel_1 = par[1],
            covariavel_2 = par[2],
            correlacao = cor(
              .x[[par[1]]],
              .x[[par[2]]],
              use = "complete.obs"
            )
          )
        }
      )
    ) |>
    ungroup()
)


# Médias mensais

series_mensais_covariaveis <-
  covariaveis_treinamento |>
  mutate(
    data_mes = as.Date(
      format(data, "%Y-%m-01")
    )
  ) |>
  group_by(
    localidade,
    clima,
    covariavel,
    data_mes
  ) |>
  summarise(
    valor_medio = mean(valor),
    .groups = "drop"
  )


# Perfil sazonal

perfil_sazonal_covariaveis <-
  covariaveis_treinamento |>
  group_by(
    localidade,
    clima,
    covariavel,
    mes
  ) |>
  summarise(
    valor_medio = mean(valor),
    .groups = "drop"
  )


# Séries mensais das localidades selecionadas

grafico_series_covariaveis_selecionadas <-
  series_mensais_covariaveis |>
  filter(
    localidade %in% localidades_selecionadas
  ) |>
  ggplot(
    aes(
      x = data_mes,
      y = valor_medio,
      color = clima
    )
  ) +
  geom_line(linewidth = 0.4) +
  facet_grid(
    covariavel ~ localidade,
    scales = "free_y",
    labeller = labeller(
      covariavel = rotulos_covariaveis
    )
  ) +
  scale_x_date(
    date_breaks = "5 years",
    date_labels = "%Y"
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = "Condição climática"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "bottom"
  )

grafico_series_covariaveis_selecionadas


# Séries mensais completas

grafico_series_covariaveis_completas <-
  ggplot(
    series_mensais_covariaveis,
    aes(
      x = data_mes,
      y = valor_medio,
      color = clima
    )
  ) +
  geom_line(linewidth = 0.35) +
  facet_grid(
    covariavel ~ localidade,
    scales = "free_y",
    labeller = labeller(
      covariavel = rotulos_covariaveis
    )
  ) +
  scale_x_date(
    date_breaks = "5 years",
    date_labels = "%Y"
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = "Condição climática"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "bottom"
  )

grafico_series_covariaveis_completas


# Perfil sazonal

grafico_perfil_sazonal_covariaveis <-
  ggplot(
    perfil_sazonal_covariaveis,
    aes(
      x = mes,
      y = valor_medio,
      color = localidade
    )
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.1) +
  facet_grid(
    covariavel ~ clima,
    scales = "free_y",
    labeller = labeller(
      covariavel = rotulos_covariaveis
    )
  ) +
  scale_x_continuous(
    breaks = 1:12
  ) +
  labs(
    x = "Mês",
    y = NULL,
    color = "Localidade"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom"
  )

grafico_perfil_sazonal_covariaveis


# Organizar correlações para o gráfico

correlacoes_grafico <-
  correlacoes_covariaveis_alvos |>
  pivot_longer(
    cols = starts_with("correlacao_"),
    names_to = "analise",
    values_to = "correlacao"
  ) |>
  mutate(
    localidade = factor(
      localidade,
      levels = ordem_localidades
    ),
    covariavel = factor(
      covariavel,
      levels = rev(covariaveis)
    ),
    alvo = factor(
      alvo,
      levels = alvos
    ),
    analise = recode(
      analise,
      correlacao_nivel = "Nível atual",
      correlacao_variacao_atual = "Variação no mesmo dia",
      correlacao_variacao_futura = "Variação no dia seguinte"
    ),
    analise = factor(
      analise,
      levels = c(
        "Nível atual",
        "Variação no mesmo dia",
        "Variação no dia seguinte"
      )
    ),
    rotulo = str_replace(
      sprintf("%.2f", correlacao),
      "\\.",
      ","
    ),
    texto_claro = abs(correlacao) > 0.55
  )


# Mapa das correlações

grafico_correlacoes_covariaveis_alvos <-
  ggplot(
    correlacoes_grafico,
    aes(
      x = localidade,
      y = covariavel,
      fill = correlacao
    )
  ) +
  geom_tile(color = "white") +
  geom_text(
    aes(
      label = rotulo,
      color = texto_claro
    ),
    size = 3
  ) +
  facet_grid(
    analise ~ alvo
  ) +
  scale_y_discrete(
    labels = rotulos_covariaveis
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  scale_color_manual(
    values = c(
      "FALSE" = "black",
      "TRUE" = "white"
    ),
    guide = "none"
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Correlação"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    legend.position = "bottom"
  )

grafico_correlacoes_covariaveis_alvos


# Pastas

dir.create(
  "tabelas",
  showWarnings = FALSE
)

dir.create(
  "figuras",
  showWarnings = FALSE
)


# Exportação das tabelas

write_csv(
  resumo_covariaveis,
  "tabelas/resumo_covariaveis.csv"
)

write_csv(
  resumo_precipitacao,
  "tabelas/resumo_precipitacao.csv"
)

write_csv(
  correlacoes_covariaveis_alvos,
  "tabelas/correlacoes_covariaveis_alvos.csv"
)

write_csv(
  correlacoes_entre_covariaveis,
  "tabelas/correlacoes_entre_covariaveis.csv"
)


# Exportação dos gráficos

ggsave(
  "figuras/series_mensais_covariaveis_selecionadas.pdf",
  grafico_series_covariaveis_selecionadas,
  width = 12,
  height = 10
)

ggsave(
  "figuras/series_mensais_covariaveis_completas.pdf",
  grafico_series_covariaveis_completas,
  width = 15,
  height = 10
)

ggsave(
  "figuras/perfil_sazonal_covariaveis.pdf",
  grafico_perfil_sazonal_covariaveis,
  width = 12,
  height = 10
)

ggsave(
  "figuras/correlacoes_covariaveis_alvos.pdf",
  grafico_correlacoes_covariaveis_alvos,
  width = 16,
  height = 18,
  units = "cm"
)



# As covariáveis apresentam padrões sazonais claros, mas com diferenças
# importantes entre as localidades

# A precipitação média e a frequência de dias com mais de 1 mm
# diminuem das localidades úmidas para as mais secas

# Valores positivos muito pequenos de precipitação são frequentes
# Por isso, o limite de 1 mm é mais informativo do que apenas separar
# dias com precipitação igual ou diferente de zero

# A umidade relativa apresenta forte correlação positiva com os níveis
# de GWETTOP e GWETROOT

# As correlações em nível podem refletir sazonalidade compartilhada
# e não demonstram ganho preditivo adicional

# A precipitação apresenta a associação mais consistente com a mudança
# diária das duas variáveis-alvo

# Para GWETTOP, a relação entre a precipitação atual e a mudança no dia
# seguinte é fraca e varia entre as localidades

# Para GWETROOT, a precipitação mantém associação positiva com a mudança
# no dia seguinte, sugerindo uma resposta mais lenta da zona radicular

# Radiação e evapotranspiração apresentam associação negativa com as
# mudanças diárias da umidade

# Algumas covariáveis são fortemente correlacionadas, principalmente
# temperatura e umidade relativa em Belém e Manaus

# A inclusão das covariáveis deverá ser avaliada fora da amostra e
# comparada com modelos que utilizam apenas o histórico da umidade
