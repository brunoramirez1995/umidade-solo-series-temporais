# Gráficos de Defasagem das Variáveis Alvo

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


# Construir as nove defasagens

dados_defasagens <- map_dfr(
  1:9,
  function(k) {
    
    alvos_treinamento |>
      group_by(localidade, clima, variavel) |>
      arrange(data, .by_group = TRUE) |>
      transmute(
        data,
        umidade_atual = umidade,
        umidade_defasada = lag(umidade, n = k),
        k = k
      ) |>
      ungroup() |>
      drop_na(umidade_atual, umidade_defasada)
  }
)


# Rótulos dos painéis

dados_defasagens <- dados_defasagens |>
  mutate(
    defasagem = factor(
      k,
      levels = 1:9,
      labels = paste("k =", 1:9)
    )
  )


head(dados_defasagens)



# Correlação da série com suas defasagens

resumo_correlacoes_defasagens <- dados_defasagens |>
  group_by(
    localidade,
    clima,
    variavel,
    k
  ) |>
  summarise(
    correlacao = cor(
      umidade_atual,
      umidade_defasada
    ),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = k,
    values_from = correlacao,
    names_prefix = "cor_lag_"
  )


resumo_correlacoes_defasagens



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


# Cores das condições climáticas

cores_clima <- c(
  "Úmido" = "#619CFF",
  "Tropical sazonal" = "#00BA38",
  "Mais seco" = "#F8766D"
)


# Função para construir cada página

construir_grafico_defasagem <- function(local, alvo) {
  
  dados_grafico <- dados_defasagens |>
    filter(
      localidade == local,
      variavel == alvo
    )
  
  clima_local <- first(dados_grafico$clima)
  cor_local <- unname(cores_clima[clima_local])
  
  ggplot(
    dados_grafico,
    aes(
      x = umidade_defasada,
      y = umidade_atual
    )
  ) +
    geom_point(
      color = cor_local,
      alpha = 0.08,
      size = 0.25
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.4
    ) +
    facet_wrap(
      ~ defasagem,
      ncol = 3
    ) +
    coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    labs(
      title = paste(local, "-", alvo),
      subtitle = paste("Condição climática:", clima_local),
      x = expression(Y[t-k]),
      y = expression(Y[t])
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9)
    )
}




# Pastas para os resultados

dir.create("figuras", showWarnings = FALSE)
dir.create("tabelas", showWarnings = FALSE)


# PDF com as 12 séries

pdf(
  "figuras/lagplots_series_diarias.pdf",
  width = 10,
  height = 10
)


for (local in localidades) {
  
  for (alvo in variaveis) {
    
    grafico <- construir_grafico_defasagem(
      local,
      alvo
    )
    
    print(grafico)
  }
}


dev.off()


# Salvar resumo numérico

write_csv(
  resumo_correlacoes_defasagens,
  "tabelas/resumo_correlacoes_defasagens.csv"
)


# Os gráficos mostram forte associação entre Y_t e seus valores defasados

# A correlação diminui gradualmente conforme a defasagem aumenta,
# mas permanece elevada até o atraso 9 em todas as séries

# A perda de dependência foi mais rápida nas localidades mais secas

# GWETROOT apresentou nuvens geralmente mais estreitas,
# indicando maior estabilidade diária que GWETTOP

# GWETTOP apresentou maior dispersão e diferenças entre faixas de umidade,
# principalmente em Brasília, Montes Claros e Petrolina

# Em algumas séries, aparecem concentrações distintas em níveis baixos
# e altos de umidade, compatíveis com a alternância sazonal

# Não foi observado um limiar único e claramente definido para todas as séries

# Os formatos encontrados justificam investigar dependência do estado,
# mas ainda não comprovam a presença de não linearidade

# Parte da estrutura também pode resultar da sazonalidade,
# da precipitação, dos limites físicos e do arredondamento dos dados


# Exportar os lag plots como imagens

dir.create(
  "figuras/lagplots",
  showWarnings = FALSE,
  recursive = TRUE
)


# Nomes seguros para os arquivos

nomes_localidades <- c(
  "Belém" = "belem",
  "Manaus" = "manaus",
  "Brasília" = "brasilia",
  "Uberlândia" = "uberlandia",
  "Montes Claros" = "montes_claros",
  "Petrolina" = "petrolina"
)


# Uma imagem para cada série

for (local in localidades) {
  
  for (alvo in variaveis) {
    
    grafico <- construir_grafico_defasagem(
      local,
      alvo
    )
    
    nome_arquivo <- paste0(
      "figuras/lagplots/lagplots_",
      nomes_localidades[[local]],
      "_",
      str_to_lower(alvo),
      ".png"
    )
    
    ggsave(
      filename = nome_arquivo,
      plot = grafico,
      width = 8,
      height = 8,
      dpi = 300,
      bg = "white"
    )
  }
}


# Lag plots selecionados para o artigo
# Brasília, duas variáveis e defasagens 1, 5 e 9

dados_lagplots_brasilia <- dados_defasagens |>
  filter(
    localidade == "Brasília",
    k %in% c(1, 5, 9)
  ) |>
  mutate(
    defasagem = factor(
      k,
      levels = c(1, 5, 9),
      labels = c("k = 1", "k = 5", "k = 9")
    )
  ) |>
  group_by(variavel, k) |>
  filter((row_number() - 1) %% 3 == 0) |>
  ungroup()

grafico_lagplots_brasilia <- ggplot(
  dados_lagplots_brasilia,
  aes(
    x = umidade_defasada,
    y = umidade_atual,
    color = variavel
  )
) +
  geom_point(
    alpha = 0.12,
    size = 0.35
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.4
  ) +
  facet_grid(
    rows = vars(variavel),
    cols = vars(defasagem)
  ) +
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  scale_color_manual(
    values = c(
      "GWETTOP" = "#D55E00",
      "GWETROOT" = "#0072B2"
    ),
    guide = "none"
  ) +
  labs(
    x = expression(Y[t-k]),
    y = expression(Y[t])
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = "figuras/lagplots_brasilia_alvos.pdf",
  plot = grafico_lagplots_brasilia,
  width = 10,
  height = 6.5,
  bg = "white"
)
