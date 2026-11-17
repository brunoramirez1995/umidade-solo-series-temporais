

# Decomposição das Variáveis Alvo

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


# Média mensal de GWETTOP

serie_mensal <- dados_treinamento |>
  mutate(
    data_mes = as.Date(format(data, "%Y-%m-01"))
  ) |>
  group_by(localidade, data_mes) |>
  summarise(
    umidade_media = mean(GWETTOP),
    .groups = "drop"
  )


# Série de Brasília

brasilia <- serie_mensal |>
  filter(localidade == "Brasília") |>
  arrange(data_mes) |>
  pull(umidade_media)


# Transformar em série temporal mensal

brasilia_ts <- ts(
  brasilia,
  start = c(2001, 1),
  frequency = 12
)

brasilia_ts

start(brasilia_ts)
end(brasilia_ts)
frequency(brasilia_ts)



# Decomposição aditiva

decomposicao_aditiva <- decompose(
  brasilia_ts,
  type = "additive"
)

plot(decomposicao_aditiva)


# Decomposição multiplicativa

decomposicao_multiplicativa <- decompose(
  brasilia_ts,
  type = "multiplicative"
)

plot(decomposicao_multiplicativa)


# Erros da decomposição aditiva

erro_aditivo <- brasilia_ts -
  decomposicao_aditiva$trend -
  decomposicao_aditiva$seasonal


# Erros da decomposição multiplicativa

erro_multiplicativo <- brasilia_ts -
  decomposicao_multiplicativa$trend *
  decomposicao_multiplicativa$seasonal


# Comparação

comparacao_decomposicao <- tibble(
  tipo = c("Aditiva", "Multiplicativa"),
  rmse = c(
    sqrt(mean(erro_aditivo^2, na.rm = TRUE)),
    sqrt(mean(erro_multiplicativo^2, na.rm = TRUE))
  )
)

comparacao_decomposicao


# A decomposição aditiva apresentou menor RMSE que a multiplicativa

# A diferença foi pequena, mas o gráfico também indica uma amplitude sazonal
# relativamente estável ao longo do período

# A série apresenta um ciclo anual bem definido

# A tendência varia ao longo do tempo e não apresenta comportamento monotônico

# A componente aleatória ainda deve ser analisada pela ACF

# Para GWETTOP de Brasília, a decomposição aditiva foi considerada mais adequada



# Comparação para todas as séries


# Médias mensais das duas variáveis

serie_mensal_alvos <- dados_treinamento |>
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
    data_mes = as.Date(format(data, "%Y-%m-01"))
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


# Comparar as decomposições

comparar_decomposicao <- function(local, alvo) {
  
  dados_serie <- serie_mensal_alvos |>
    filter(
      localidade == local,
      variavel == alvo
    ) |>
    arrange(data_mes)
  
  serie_ts <- ts(
    dados_serie$umidade_media,
    start = c(2001, 1),
    frequency = 12
  )
  
  aditiva <- decompose(
    serie_ts,
    type = "additive"
  )
  
  multiplicativa <- decompose(
    serie_ts,
    type = "multiplicative"
  )
  
  erro_aditivo <- serie_ts -
    aditiva$trend -
    aditiva$seasonal
  
  erro_multiplicativo <- serie_ts -
    multiplicativa$trend *
    multiplicativa$seasonal
  
  rmse_aditivo <- sqrt(
    mean(erro_aditivo^2, na.rm = TRUE)
  )
  
  rmse_multiplicativo <- sqrt(
    mean(erro_multiplicativo^2, na.rm = TRUE)
  )
  
  tibble(
    localidade = local,
    clima = first(dados_serie$clima),
    variavel = alvo,
    rmse_aditivo = rmse_aditivo,
    rmse_multiplicativo = rmse_multiplicativo,
    diferenca_rmse = rmse_multiplicativo - rmse_aditivo
  )
}


# Combinações de localidade e variável

combinacoes <- serie_mensal_alvos |>
  distinct(localidade, variavel)


# Resultado das 12 séries

(resultado_decomposicoes <- map2_dfr(
  combinacoes$localidade,
  combinacoes$variavel,
  comparar_decomposicao
) |>
    arrange(clima, localidade, variavel))


# Diferença positiva: a decomposição aditiva apresentou menor RMSE

# Diferença negativa: a decomposição multiplicativa apresentou menor RMSE

# Valores próximos de zero indicam pouca diferença entre os métodos






# Para GWETTOP, a decomposição aditiva apresentou menor RMSE
# em cinco das seis localidades

# Petrolina foi a única localidade em que a decomposição
# multiplicativa apresentou menor RMSE para GWETTOP

# Para GWETROOT, a decomposição multiplicativa apresentou
# menor RMSE em cinco das seis localidades

# Belém foi a única localidade em que a decomposição
# aditiva apresentou menor RMSE para GWETROOT

# As diferenças entre os métodos foram pequenas na maioria das séries

# A escolha final também deve considerar o comportamento visual
# das componentes e não apenas o RMSE

# Essa comparação avalia a representação da série e não
# o desempenho de modelos de previsão




# Gráficos das decomposições


# Ordem dos gráficos

localidades_decomposicao <- c(
  "Belém",
  "Manaus",
  "Brasília",
  "Uberlândia",
  "Montes Claros",
  "Petrolina"
)

variaveis_decomposicao <- c(
  "GWETTOP",
  "GWETROOT"
)

tipos_decomposicao <- c(
  "additive",
  "multiplicative"
)


# Pasta dos gráficos

dir.create(
  "figuras",
  showWarnings = FALSE
)


# PDF com todas as comparações

pdf(
  "figuras/comparacao_grafica_decomposicoes.pdf",
  width = 10,
  height = 7
)

for (local in localidades_decomposicao) {
  
  for (alvo in variaveis_decomposicao) {
    
    dados_serie <- serie_mensal_alvos |>
      filter(
        localidade == local,
        variavel == alvo
      ) |>
      arrange(data_mes)
    
    serie_ts <- ts(
      dados_serie$umidade_media,
      start = c(2001, 1),
      frequency = 12
    )
    
    for (tipo in tipos_decomposicao) {
      
      decomposicao <- decompose(
        serie_ts,
        type = tipo
      )
      
      plot(
        decomposicao,
        xlab = paste(local, "-", alvo)
      )
    }
  }
}

dev.off()

# Foram geradas as decomposições aditiva e multiplicativa
# para as 12 séries mensais

# O arquivo possui 24 páginas e será utilizado para comparar
# visualmente a estabilidade das componentes

# Os gráficos finais serão definidos após a comparação
# entre o RMSE e o comportamento das decomposições



# Os gráficos confirmam a presença de sazonalidade anual em todas as séries

# A tendência estimada é semelhante nas decomposições aditiva e multiplicativa

# As principais diferenças aparecem na escala das componentes
# sazonal e aleatória

# GWETTOP apresentou comportamento predominantemente aditivo

# Petrolina foi a única localidade com indicação multiplicativa
# para GWETTOP

# GWETROOT apresentou comportamento predominantemente multiplicativo

# Belém foi a única localidade com indicação aditiva para GWETROOT

# As diferenças visuais entre os métodos foram pequenas,
# em concordância com os valores próximos de RMSE

# A componente aleatória ainda não pode ser considerada ruído branco

# A dependência temporal restante será investigada pela ACF





# Tipo adotado para cada variável

resultado_decomposicoes <- resultado_decomposicoes |>
  mutate(
    tipo_adotado = recode(
      variavel,
      GWETTOP = "Aditiva",
      GWETROOT = "Multiplicativa"
    )
  )

resultado_decomposicoes

# A escolha não foi realizada separadamente por localidade

# Foi adotado o mesmo tipo de decomposição dentro de cada variável
# para manter a comparação entre os domínios climáticos

# A decomposição aditiva foi adotada para GWETTOP

# A decomposição multiplicativa foi adotada para GWETROOT

# Cada escolha apresentou menor RMSE em cinco das seis localidades

# Petrolina e Belém foram exceções, mas as diferenças de RMSE
# não justificaram alterar o método apenas nessas localidades

write_csv(
  resultado_decomposicoes,
  "tabelas/comparacao_decomposicoes.csv")


# Figura para o artigo
# Figuras principais do artigo


# GWETTOP de Brasília

dados_brasilia_gwettop <- serie_mensal_alvos |>
  filter(
    localidade == "Brasília",
    variavel == "GWETTOP"
  ) |>
  arrange(data_mes)


componentes_brasilia_gwettop <- tibble(
  data = dados_brasilia_gwettop$data_mes,
  observada = as.numeric(brasilia_ts),
  tendencia = as.numeric(decomposicao_aditiva$trend),
  sazonalidade = as.numeric(decomposicao_aditiva$seasonal),
  irregular = as.numeric(decomposicao_aditiva$random)
) |>
  pivot_longer(
    cols = -data,
    names_to = "componente",
    values_to = "valor"
  ) |>
  mutate(
    componente = factor(
      componente,
      levels = c(
        "observada",
        "tendencia",
        "sazonalidade",
        "irregular"
      ),
      labels = c(
        "Série observada",
        "Tendência",
        "Componente sazonal",
        "Componente irregular"
      )
    )
  ) |>
  filter(!is.na(valor))


grafico_decomposicao_gwettop <- ggplot(
  componentes_brasilia_gwettop,
  aes(x = data, y = valor)
) +
  geom_line(
    linewidth = 0.4,
    color = "#D55E00"
  ) +
  facet_wrap(
    ~ componente,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_x_date(
    date_breaks = "4 years",
    date_labels = "%Y"
  ) +
  labs(
    x = "Ano",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

grafico_decomposicao_gwettop


ggsave(
  "figuras/decomposicao_brasilia_gwettop.pdf",
  grafico_decomposicao_gwettop,
  width = 10,
  height = 7
)



# GWETROOT de Brasília

dados_brasilia_gwetroot <- serie_mensal_alvos |>
  filter(
    localidade == "Brasília",
    variavel == "GWETROOT"
  ) |>
  arrange(data_mes)


brasilia_gwetroot_ts <- ts(
  dados_brasilia_gwetroot$umidade_media,
  start = c(2001, 1),
  frequency = 12
)


decomposicao_brasilia_gwetroot <- decompose(
  brasilia_gwetroot_ts,
  type = "multiplicative"
)


componentes_brasilia_gwetroot <- tibble(
  data = dados_brasilia_gwetroot$data_mes,
  observada = as.numeric(brasilia_gwetroot_ts),
  tendencia = as.numeric(decomposicao_brasilia_gwetroot$trend),
  sazonalidade = as.numeric(decomposicao_brasilia_gwetroot$seasonal),
  irregular = as.numeric(decomposicao_brasilia_gwetroot$random)
) |>
  pivot_longer(
    cols = -data,
    names_to = "componente",
    values_to = "valor"
  ) |>
  mutate(
    componente = factor(
      componente,
      levels = c(
        "observada",
        "tendencia",
        "sazonalidade",
        "irregular"
      ),
      labels = c(
        "Série observada",
        "Tendência",
        "Fator sazonal",
        "Fator irregular"
      )
    )
  ) |>
  filter(!is.na(valor))


grafico_decomposicao_gwetroot <- ggplot(
  componentes_brasilia_gwetroot,
  aes(x = data, y = valor)
) +
  geom_line(
    linewidth = 0.4,
    color = "#0072B2"
  ) +
  facet_wrap(
    ~ componente,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_x_date(
    date_breaks = "4 years",
    date_labels = "%Y"
  ) +
  labs(
    x = "Ano",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

grafico_decomposicao_gwetroot


ggsave(
  "figuras/decomposicao_brasilia_gwetroot.pdf",
  grafico_decomposicao_gwetroot,
  width = 10,
  height = 7
)

# Brasília foi utilizada para apresentar as duas decomposições
# no texto principal

# GWETTOP foi decomposta pela estrutura aditiva

# GWETROOT foi decomposta pela estrutura multiplicativa

# Na decomposição aditiva, as componentes sazonal e irregular
# são centradas em zero

# Na decomposição multiplicativa, os fatores sazonal e irregular
# são centrados em um

# As demais localidades serão apresentadas no apêndice