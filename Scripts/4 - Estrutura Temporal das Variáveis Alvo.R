

# Estrutura Temporal das Variáveis Alvo

library(tidyverse)
library(readr)
library(tidyr)
library(stringr)
library(purrr)


# Ler arquivo

dados <- read_csv(
  "Dados/nasa_power_diario.csv",
  show_col_types = FALSE
)

dados$data <- as.Date(dados$data)


# Conjunto de treinamento

dados_treinamento <- dados |>
  filter(data <= as.Date("2021-12-31"))

# O período de 2001 a 2021 será utilizado para análise e ajuste dos modelos

# Os anos de 2022 a 2024 serão reservados para validação

# O ano de 2025 será utilizado como teste final


# Organizar as variáveis alvo

alvos_treinamento <- dados_treinamento |>
  select(
    localidade,
    clima,
    data,
    ano,
    mes,
    GWETTOP,
    GWETROOT
  ) |>
  pivot_longer(
    cols = c(GWETTOP, GWETROOT),
    names_to = "variavel",
    values_to = "umidade"
  )


# Ordem das localidades

alvos_treinamento$localidade <- factor(
  alvos_treinamento$localidade,
  levels = c(
    "Belém",
    "Manaus",
    "Brasília",
    "Uberlândia",
    "Montes Claros",
    "Petrolina"
  )
)

alvos_treinamento$clima <- factor(
  alvos_treinamento$clima,
  levels = c(
    "Úmido",
    "Tropical sazonal",
    "Mais seco"
  )
)

alvos_treinamento$variavel <- factor(
  alvos_treinamento$variavel,
  levels = c(
    "GWETTOP",
    "GWETROOT"
  )
)






# As médias mensais facilitam a visualização da sazonalidade

# As médias anuais ajudam a observar mudanças de longo prazo

# Essas agregações são apenas descritivas; a série original permanece diária


# Série mensal

serie_mensal_treinamento <- alvos_treinamento |>
  mutate(
    data_mes = as.Date(
      format(data, "%Y-%m-01")
    )
  ) |>
  group_by(
    localidade,
    clima,
    variavel,
    data_mes
  ) |>
  summarise(
    umidade_media = mean(umidade),
    .groups = "drop"
  )


# Série anual

serie_anual_treinamento <- alvos_treinamento |>
  group_by(
    localidade,
    clima,
    variavel,
    ano
  ) |>
  summarise(
    umidade_media = mean(umidade),
    .groups = "drop"
  )


# Resumo da tendência

(resumo_tendencia <- serie_anual_treinamento |>
    group_by(
      localidade,
      clima,
      variavel
    ) |>
    summarise(
      media_2001 = umidade_media[ano == 2001],
      media_2021 = umidade_media[ano == 2021],
      diferenca_periodo = media_2021 - media_2001,
      inclinacao_decada =
        unname(
          coef(
            lm(umidade_media ~ ano)
          )[2]
        ) * 10,
      .groups = "drop"
    )
)

# Belém e Manaus apresentam inclinação positiva no período

# Brasília e Uberlândia apresentam redução nas duas camadas

# Nas localidades mais secas, a redução é mais clara em GWETTOP

# GWETROOT apresenta maior estabilidade em Montes Claros e Petrolina

# A inclinação resume a direção geral, mas não confirma uma tendência estatística


# Perfil médio mensal

perfil_sazonal_treinamento <- alvos_treinamento |>
  group_by(
    localidade,
    clima,
    variavel,
    mes
  ) |>
  summarise(
    umidade_media = mean(umidade),
    .groups = "drop"
  )


# Resumo da sazonalidade

(resumo_sazonal <- perfil_sazonal_treinamento |>
    group_by(
      localidade,
      clima,
      variavel
    ) |>
    summarise(
      menor_media = min(umidade_media),
      mes_menor_media = mes[
        which.min(umidade_media)
      ],
      maior_media = max(umidade_media),
      mes_maior_media = mes[
        which.max(umidade_media)
      ],
      amplitude_sazonal =
        maior_media - menor_media,
      .groups = "drop"
    )
)


# Os menores níveis de umidade geralmente ocorrem entre agosto e outubro

# O mês de maior umidade varia entre as localidades

# GWETTOP apresenta maior amplitude sazonal nas localidades sazonais e mais secas

# Brasília apresenta a maior amplitude sazonal de GWETTOP

# GWETROOT apresenta pouca variação sazonal em Petrolina

# Localidades do mesmo grupo climático não possuem comportamento idêntico...


# Gráfico de tendência anual


grafico_tendencia_anual <- ggplot(
  serie_anual_treinamento,
  aes(
    x = ano,
    y = umidade_media,
    color = clima
  )
) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  facet_grid(
    variavel ~ localidade
  ) +
  scale_x_continuous(
    breaks = seq(2001, 2021, 5)
  ) +
  labs(
    x = "Ano",
    y = "Umidade média anual",
    color = "Condição climática"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "bottom"
  )

grafico_tendencia_anual

# Localidades

localidades_selecionadas <- c(
  "Belém",
  "Brasília",
  "Petrolina"
)

serie_mensal_selecionadas <- serie_mensal_treinamento |>
  filter(
    localidade %in% localidades_selecionadas
  )


grafico_mensal_selecionadas <- ggplot(
  serie_mensal_selecionadas,
  aes(
    x = data_mes,
    y = umidade_media,
    color = variavel
  )
) +
  geom_line(linewidth = 0.5) +
  facet_wrap(
    ~ localidade,
    nrow = 1
  ) +
  scale_x_date(
    date_breaks = "4 years",
    date_labels = "%Y"
  ) +
  labs(
    x = NULL,
    y = "Umidade média mensal",
    color = "Variável"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "bottom"
  )

grafico_mensal_selecionadas


grafico_sazonal_treinamento <- ggplot(
  perfil_sazonal_treinamento,
  aes(
    x = mes,
    y = umidade_media,
    color = localidade
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.3) +
  facet_grid(
    variavel ~ clima
  ) +
  scale_x_continuous(
    breaks = 1:12
  ) +
  labs(
    x = "Mês",
    y = "Umidade média",
    color = "Localidade"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom"
  )

grafico_sazonal_treinamento


# Pastas

dir.create(
  "tabelas",
  showWarnings = FALSE
)

dir.create(
  "figuras",
  showWarnings = FALSE
)


# Tabelas

write_csv(
  resumo_tendencia,
  "tabelas/resumo_tendencia.csv"
)

write_csv(
  resumo_sazonal,
  "tabelas/resumo_sazonal.csv"
)


# Gráficos

ggsave(
  "figuras/tendencia_anual.pdf",
  grafico_tendencia_anual,
  width = 13,
  height = 6
)

ggsave(
  "figuras/serie_mensal_selecionadas.pdf",
  grafico_mensal_selecionadas,
  width = 11,
  height = 4.5
)

ggsave(
  "figuras/perfil_sazonal_treinamento.pdf",
  grafico_sazonal_treinamento,
  width = 12,
  height = 7
)



# Tendência anual das localidades selecionadas

tendencia_anual_selecionadas <- serie_anual_treinamento |>
  filter(
    localidade %in% localidades_selecionadas
  )


grafico_tendencia_anual_selecionadas <- ggplot(
  tendencia_anual_selecionadas,
  aes(
    x = ano,
    y = umidade_media,
    color = variavel
  )
) +
  geom_line(
    linewidth = 0.6,
    alpha = 0.7
  ) +
  geom_point(
    size = 1.2,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  facet_wrap(
    ~ localidade,
    nrow = 1
  ) +
  scale_x_continuous(
    breaks = seq(2001, 2021, 5)
  ) +
  labs(
    x = "Ano",
    y = "Umidade média anual",
    color = "Variável"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "bottom"
  )

grafico_tendencia_anual_selecionadas

ggsave(
  "figuras/tendencia_anual_selecionadas.pdf",
  grafico_tendencia_anual_selecionadas,
  width = 11,
  height = 4.5
)
