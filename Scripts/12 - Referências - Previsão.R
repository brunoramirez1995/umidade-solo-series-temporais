
# Persistência e previsão sazonal ingênua: um dia à frente

# Pacotes e configurações

library(tidyverse)

source("Scripts/apoio_modelagem.R", encoding = "UTF-8")

arquivo_dados <- "Dados/nasa_power_diario.csv"
pasta_resultados <- "resultados/modelagem/referencias"


# Leitura das colunas necessárias

dados <- read_csv(
  arquivo_dados,
  col_types = cols_only(
    localidade_id = col_character(),
    localidade = col_character(),
    clima = col_character(),
    data = col_date(),
    GWETTOP = col_double(),
    GWETROOT = col_double()
  ),
  show_col_types = FALSE
)


# 2025 fica fora desta avaliação

dados <- dados |>
  filter(data >= inicio_treinamento, data <= fim_validacao) |>
  arrange(localidade_id, data)



# Previsões: cada linha corresponde à data que será prevista

# Persistência: valor observado no dia anterior.
# Sazonal: valor observado 365 dias antes da data prevista.
# Calcular os lags antes de selecionar a validação preserva o histórico.

previsoes <- dados |>
  pivot_longer(
    all_of(alvos),
    names_to = "variavel",
    values_to = "observado"
  ) |>
  arrange(localidade_id, variavel, data) |>
  group_by(localidade_id, variavel) |>
  mutate(
    data_origem = data - 1,
    persistencia = dplyr::lag(observado, 1),
    sazonal_365 = dplyr::lag(observado, 365)
  ) |>
  ungroup() |>
  filter(data >= inicio_validacao, data <= fim_validacao) |>
  rename(data_previsao = data) |>
  pivot_longer(
    c(persistencia, sazonal_365),
    names_to = "modelo",
    values_to = "previsto"
  ) |>
  mutate(
    ano = as.integer(format(data_previsao, "%Y")),
    erro = observado - previsto
  )

# A origem avança diariamente.
# Observações anteriores da validação já estão disponíveis quando a previsão do próximo dia é produzida.




# Métricas por série e por ano

metricas <- previsoes |>
  group_by(localidade_id, localidade, clima, variavel, modelo) |>
  summarise(
    calcular_metricas(observado, previsto),
    .groups = "drop"
  )

metricas_anuais <- previsoes |>
  group_by(localidade_id, localidade, clima, variavel, modelo, ano) |>
  summarise(
    calcular_metricas(observado, previsto),
    .groups = "drop"
  )


# RMSE e MAE menores indicam previsões melhores.
# O resumo total usa todos os erros, não a média dos RMSE anuais.

print(calendario_validacao, width = Inf)
print(metricas, n = Inf, width = Inf)
print(metricas_anuais, n = Inf, width = Inf)

# Interpretação dos resultados — validação de 2022 a 2024


# A persistência apresentou RMSE e MAE menores nas 12 séries.
# A vantagem também ocorreu nas 36 comparações por série e ano.

# Para prever um dia à frente, o valor do dia anterior foi mais
# informativo do que o valor observado 365 dias antes.


# Os erros da persistência foram menores em GWETROOT, acompanhando
# sua menor variação diária. A precisão dos registros também deve
# ser considerada ao interpretar os dias com erro igual a zero.


# O desempenho da referência sazonal não descarta a utilidade de
# termos harmônicos nos próximos modelos.


# O viés próximo de zero pode esconder compensação entre erros.
# A comparação deve considerar principalmente RMSE e MAE.


# A persistência será a principal referência para avaliar o ganho
# dos próximos modelos. Estes resultados ainda não incluem 2025.


# Comparação visual

# Cada camada tem sua própria escala de RMSE.

grafico <- ggplot(
  metricas,
  aes(x = RMSE, y = localidade, color = modelo)
) +
  geom_point(size = 2.5) +
  facet_wrap(~ variavel, scales = "free_x") +
  scale_color_manual(
    values = c(
      persistencia = "#2166AC",
      sazonal_365 = "#B35806"
    ),
    breaks = c("persistencia", "sazonal_365"),
    labels = c("Persistência", "Sazonal (365 dias)")
  ) +
  labs(
    x = "RMSE na validação (2022-2024)",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

if (interactive()) print(grafico)


# Salvar resultados para comparar com as próximas famílias

dir.create(
  pasta_resultados,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figuras/principais",
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  previsoes,
  file.path(pasta_resultados, "previsoes_validacao.csv")
)

write_csv(
  metricas,
  file.path(pasta_resultados, "metricas_validacao.csv")
)

write_csv(
  metricas_anuais,
  file.path(pasta_resultados, "metricas_validacao_anuais.csv")
)

ggsave(
  "figuras/principais/referencias_validacao.pdf",
  grafico,
  width = 9,
  height = 4.5
)

