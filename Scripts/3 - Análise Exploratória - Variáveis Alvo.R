
# Análise Exploratória

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


# Divisão temporal

fim_treinamento <- as.Date("2021-12-31")
fim_validacao <- as.Date("2024-12-31")

dados <- dados |>
  mutate(
    periodo = case_when(
      data <= fim_treinamento ~ "Treinamento",
      data <= fim_validacao ~ "Validação",
      data > fim_validacao ~ "Teste"
    ),
    periodo = factor(
      periodo,
      levels = c("Treinamento", "Validação", "Teste")
    )
  )


# Quantidade de observações em cada período

(resumo_periodos <- dados |>
    count(periodo) |>
    mutate(
      percentual = n / sum(n) * 100
    )
)



dados_treinamento <- dados |>
  filter(periodo == "Treinamento")

dados_validacao <- dados |>
  filter(periodo == "Validação")

dados_teste <- dados |>
  filter(periodo == "Teste")


# Avaliação da transformação

ajuste <- 0.005

transformacao_treinamento <- dados_treinamento |>
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
  mutate(
    umidade_ajustada = pmin(
      pmax(umidade, ajuste),
      1 - ajuste
    ),
    umidade_logit = qlogis(umidade_ajustada)
  )

# Resumo das duas escalas

(resumo_transformacao <- transformacao_treinamento |>
    group_by(localidade, clima, variavel) |>
    summarise(
      observacoes = n(),
      valores_ajustados = sum(
        umidade != umidade_ajustada
      ),
      media_original = mean(umidade),
      desvio_original = sd(umidade),
      assimetria_original =
        mean((umidade - mean(umidade))^3) /
        sd(umidade)^3,
      media_logit = mean(umidade_logit),
      desvio_logit = sd(umidade_logit),
      assimetria_logit =
        mean((umidade_logit - mean(umidade_logit))^3) /
        sd(umidade_logit)^3,
      .groups = "drop"
    )
)

write_csv(
  resumo_transformacao,
  "tabelas/resumo_transformacao.csv"
)


# Verificação da base

summary(dados)


variaveis_numericas <- c(
  "GWETTOP",
  "GWETROOT",
  "PRECTOTCORR",
  "T2M",
  "RH2M",
  "ALLSKY_SFC_SW_DWN",
  "EVPTRNS"
)




qualidade_base <- tibble(
  observacoes = nrow(dados),
  localidades = n_distinct(dados$localidade_id),
  data_inicial = min(dados$data),
  data_final = max(dados$data),
  duplicata = sum(duplicated(dados[c("localidade_id", "data")])),
  valores_ausentes = sum(is.na(dados)),
  
  codigos_menos_999 = sum(
    dados[variaveis_numericas] == -999,
    na.rm = TRUE
  ),
  
  gwettop_fora_intervalo = sum(
    dados$GWETTOP < 0 | dados$GWETTOP > 1,
    na.rm = TRUE
  ),
  
  gwetroot_fora_intervalo = sum(
    dados$GWETROOT < 0 | dados$GWETROOT > 1,
    na.rm = TRUE
  )
)

qualidade_base




# Organizando variáveis alvo

alvos <- dados_treinamento |>
  
  select(
    localidade_id,
    localidade,
    clima,
    data,
    mes,
    GWETTOP,
    GWETROOT
  ) |>
  
  pivot_longer(
    cols = c(GWETTOP, GWETROOT),
    names_to = "variavel",
    values_to = "umidade"
  )


head(alvos)

# Resumo por localidade

resumo_localidades <- alvos |>
  
  group_by(localidade, clima, variavel) |>
  
  summarise(
    observacoes = n(),
    media = mean(umidade),
    desvio_padrao = sd(umidade),
    minimo = min(umidade),
    mediana = median(umidade),
    maximo = max(umidade),
    .groups = "drop"
  )

resumo_localidades


# Todas as localidades possuem 7670 observações

# GWETTOP apresenta maior variabilidade na maioria das localidades

# GWETROOT apresenta menor variação em Petrolina (DP = 0,0367)
# e Montes Claros (DP = 0,0641)

# Localidades do mesmo grupo climático parecem ter uma diferença considerável



# Resumo por condição climática

resumo_climas <- alvos |>
  
  group_by(clima, variavel) |>
  
  summarise(
    observacoes = n(),
    media = mean(umidade),
    desvio_padrao = sd(umidade),
    mediana = median(umidade),
    .groups = "drop"
  )

resumo_climas

# A umidade é maior nas localidades úmidas (ufa)

# GWETTOP diferencia melhor as condições climáticas

# GWETROOT: as médias dos grupos mais seco e tropical sazonal são próximas

# Cada grupo reúne apenas duas localidades e deve ser interpretado apenas de forma descritiva


# Esses resultados ainda não permitem concluir sobre linearidade ou não linearidade


# Correlação de GWETTOP e GWETROOT


correlacao_alvos <- dados_treinamento |>
  
  group_by(localidade, clima) |>
  summarise(
    correlacao = cor(GWETTOP, GWETROOT),
    .groups = "drop"
  )


correlacao_alvos


# GWETTOP e GWETROOT apresentam correlação positiva forte em todas as localidades

# A associação é maior em Belém e Manaus e menor em Petrolina


# Correlação elevada não significa que as variáveis tenham a mesma dinâmica


# Identificação do mês (Agrupando)

alvos$data_mes <- as.Date(format(alvos$data, "%Y-%m-01"))

# Média mensal:

serie_mensal <- alvos |>
  group_by(localidade, clima, variavel, data_mes) |>
  summarise(
    umidade_media = mean(umidade),
    .groups = "drop"
  )

serie_mensal

# Cada linha representa um mês, uma localidade e uma variável -> fazer gráfico


grafico_serie <- ggplot(
  serie_mensal,
  aes(x = data_mes, y = umidade_media, color = clima)
) +
  geom_line(linewidth = 0.4) +
  facet_grid(variavel ~ localidade) +
  scale_x_date(
    date_breaks = "5 years",
    date_labels = "%Y"
  ) +
  labs(
    x = NULL,
    y = "Umidade média mensal",
    color = "Condição climática"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

grafico_serie

# Ordenando por condição climática


# Ordem das localidades e variáveis no gráfico

serie_mensal$localidade <- factor(
  serie_mensal$localidade,
  levels = c(
    "Belém",
    "Manaus",
    "Brasília",
    "Uberlândia",
    "Montes Claros",
    "Petrolina"
  )
)

serie_mensal$variavel <- factor(
  serie_mensal$variavel,
  levels = c("GWETTOP", "GWETROOT")
)

grafico_serie <- ggplot(
  serie_mensal,
  aes(x = data_mes, y = umidade_media, color = clima)
) +
  geom_line(linewidth = 0.4) +
  facet_grid(variavel ~ localidade) +
  scale_x_date(
    date_breaks = "5 years",
    date_labels = "%Y"
  ) +
  labs(
    x = NULL,
    y = "Umidade média mensal",
    color = "Condição climática"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  )

grafico_serie


# As séries apresentam padrões sazonais ao longo do período

# GWETTOP varia mais do que GWETROOT na maioria das localidades

# GWETROOT apresenta comportamento bastante estável em Petrolina

# Localidades do mesmo grupo climático possuem dinâmicas diferentes

# Uberlândia apresenta possível redução da umidade ao longo do tempo




# Perfil médio de cada mês do ano

perfil_sazonal <- alvos |>
  group_by(localidade, clima, variavel, mes) |>
  summarise(
    umidade_media = mean(umidade),
    .groups = "drop"
  )

perfil_sazonal

perfil_sazonal$clima <- factor(
  perfil_sazonal$clima,
  levels = c("Úmido", "Tropical sazonal", "Mais seco")
)

perfil_sazonal$variavel <- factor(
  perfil_sazonal$variavel,
  levels = c("GWETTOP", "GWETROOT")
)

perfil_sazonal$localidade <- factor(
  perfil_sazonal$localidade,
  levels = c(
    "Belém",
    "Manaus",
    "Brasília",
    "Uberlândia",
    "Montes Claros",
    "Petrolina"
  )
)



# Perfil sazonal das variáveis

grafico_sazonal <- ggplot(
  perfil_sazonal,
  aes(x = mes, y = umidade_media, color = localidade)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.3) +
  facet_grid(variavel ~ clima) +
  scale_x_continuous(breaks = 1:12) +
  labs(
    x = "Mês",
    y = "Umidade média",
    color = "Localidade"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

grafico_sazonal

# As duas variáveis apresentam comportamento sazonal

# A umidade geralmente diminui no meio do ano e aumenta no final do ano

# GWETTOP apresenta maior amplitude sazonal na maioria das localidades

# Nas localidades úmidas, as duas camadas possuem perfis semelhantes

# Nas localidades mais secas, GWETROOT apresenta menor variação

# GWETROOT é muito estável em Petrolina

# Localidades do mesmo grupo climático apresentam níveis e amplitudes diferentes




# Análise da distribuição



# Ordem das localidades e variáveis

alvos$localidade <- factor(
  alvos$localidade,
  levels = c(
    "Belém",
    "Manaus",
    "Brasília",
    "Uberlândia",
    "Montes Claros",
    "Petrolina"
  )
)

alvos$variavel <- factor(
  alvos$variavel,
  levels = c("GWETTOP", "GWETROOT")
)


grafico_distribuicao <- ggplot(
  alvos,
  aes(x = localidade, y = umidade, fill = clima)
) +
  geom_boxplot(outlier.alpha = 0.08) +
  facet_wrap(~ variavel, ncol = 1) +
  labs(
    x = NULL,
    y = "Umidade do solo",
    fill = "Condição climática"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )

grafico_distribuicao

# GWETTOP apresenta valores maiores nas localidades úmidas

# A umidade superficial diminui em direção às localidades mais secas

# GWETROOT apresenta menor dispersão em Montes Claros e Petrolina

# Petrolina possui uma distribuição de GWETROOT bastante concentrada

# Os pontos externos podem representar períodos úmidos e não necessariamente outlier


# Vou escolher um ano e uma localidade X para avaliar série por dia

# Série diária de Brasília em 2021

serie_diaria_brasilia <- alvos |>
  filter(
    localidade == "Brasília",
    data >= as.Date("2021-01-01"),
    data <= as.Date("2021-12-31")
  )

grafico_diario_brasilia <- ggplot(
  serie_diaria_brasilia,
  aes(x = data, y = umidade, color = variavel)
) +
  geom_line(linewidth = 0.5) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b"
  ) +
  labs(
    title = "Umidade diária do solo em Brasília - 2021",
    x = NULL,
    y = "Umidade do solo",
    color = "Variável"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

grafico_diario_brasilia


# As duas camadas apresentam comportamento semelhante no período úmido

# Durante a estação seca, GWETTOP diminui mais rapidamente

# GWETROOT permanece mais estável entre maio e setembro

# GWETTOP apresenta recuperação acentuada no início do período úmido

# O padrão sugere possíveis regimes seco e úmido



# Salvar resultados




# Pastas para os resultados

dir.create("tabelas", showWarnings = FALSE)
dir.create("figuras", showWarnings = FALSE)


# Exportação das tabelas

write_csv(
  resumo_localidades,
  "tabelas/resumo_localidades.csv"
)

write_csv(
  resumo_climas,
  "tabelas/resumo_climas.csv"
)

write_csv(
  correlacao_alvos,
  "tabelas/correlacao_alvos.csv"
)


# Exportação dos gráficos

ggsave(
  "figuras/serie_mensal.pdf",
  grafico_serie,
  width = 13,
  height = 6
)

ggsave(
  "figuras/perfil_sazonal.pdf",
  grafico_sazonal,
  width = 12,
  height = 7
)

ggsave(
  "figuras/distribuicao_alvos.pdf",
  grafico_distribuicao,
  width = 10,
  height = 7
)

ggsave(
  "figuras/serie_diaria_brasilia_2021.pdf",
  grafico_diario_brasilia,
  width = 10,
  height = 5
)
