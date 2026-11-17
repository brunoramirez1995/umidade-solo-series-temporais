# Resposta da Umidade à Precipitação

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


# Utilizar somente o conjunto de treinamento

dados_treinamento <- dados |>
  filter(data <= as.Date("2021-12-31"))


# Organizar as duas variáveis-alvo

alvos_precipitacao <- dados_treinamento |>
  select(
    localidade,
    clima,
    data,
    PRECTOTCORR,
    GWETTOP,
    GWETROOT
  ) |>
  pivot_longer(
    cols = c(GWETTOP, GWETROOT),
    names_to = "variavel",
    values_to = "umidade"
  ) |>
  arrange(
    localidade,
    variavel,
    data
  )


# Limites utilizados para classificar o estado de umidade

limites_estado <- alvos_precipitacao |>
  group_by(
    localidade,
    clima,
    variavel
  ) |>
  summarise(
    limite_baixo = quantile(
      umidade,
      probs = 1 / 3
    ),
    limite_alto = quantile(
      umidade,
      probs = 2 / 3
    ),
    .groups = "drop"
  )

limites_estado


# Calcular as variações da umidade

alvos_precipitacao <- alvos_precipitacao |>
  left_join(
    limites_estado,
    by = c(
      "localidade",
      "clima",
      "variavel"
    )
  ) |>
  group_by(
    localidade,
    variavel
  ) |>
  mutate(
    umidade_anterior = lag(umidade),
    
    # Mudança observada no próprio dia da precipitação
    variacao_mesmo_dia =
      umidade - umidade_anterior,
    
    # Mudança observada do dia da precipitação para o dia seguinte
    variacao_dia_seguinte =
      lead(umidade) - umidade
  ) |>
  ungroup()


# Classificar o estado anterior da umidade e a precipitação

alvos_precipitacao <- alvos_precipitacao |>
  mutate(
    estado_umidade = case_when(
      umidade_anterior <= limite_baixo ~ "Baixa",
      umidade_anterior <= limite_alto ~ "Intermediária",
      umidade_anterior > limite_alto ~ "Alta",
      TRUE ~ NA_character_
    ),
    
    classe_precipitacao = case_when(
      PRECTOTCORR <= 1 ~ "Até 1 mm",
      PRECTOTCORR <= 5 ~ "1 a 5 mm",
      PRECTOTCORR > 5 ~ "Acima de 5 mm",
      TRUE ~ NA_character_
    )
  )


# Colocar os dois momentos em uma única coluna

respostas_precipitacao <- alvos_precipitacao |>
  pivot_longer(
    cols = c(
      variacao_mesmo_dia,
      variacao_dia_seguinte
    ),
    names_to = "momento",
    values_to = "variacao_umidade"
  ) |>
  mutate(
    momento = recode(
      momento,
      variacao_mesmo_dia = "Mesmo dia",
      variacao_dia_seguinte = "Dia seguinte"
    )
  ) |>
  filter(
    !is.na(estado_umidade),
    !is.na(variacao_umidade)
  )


# Definir a ordem das categorias

respostas_precipitacao <- respostas_precipitacao |>
  mutate(
    localidade = factor(
      localidade,
      levels = c(
        "Belém",
        "Manaus",
        "Brasília",
        "Uberlândia",
        "Montes Claros",
        "Petrolina"
      )
    ),
    variavel = factor(
      variavel,
      levels = c(
        "GWETTOP",
        "GWETROOT"
      )
    ),
    estado_umidade = factor(
      estado_umidade,
      levels = c(
        "Baixa",
        "Intermediária",
        "Alta"
      )
    ),
    classe_precipitacao = factor(
      classe_precipitacao,
      levels = c(
        "Até 1 mm",
        "1 a 5 mm",
        "Acima de 5 mm"
      )
    ),
    momento = factor(
      momento,
      levels = c(
        "Mesmo dia",
        "Dia seguinte"
      )
    )
  )


# Resumir a resposta à precipitação

(resumo_resposta_precipitacao <- respostas_precipitacao |>
    group_by(
      localidade,
      clima,
      variavel,
      momento,
      estado_umidade,
      classe_precipitacao
    ) |>
    summarise(
      observacoes = n(),
      variacao_media = mean(variacao_umidade),
      variacao_mediana = median(variacao_umidade),
      desvio_padrao = sd(variacao_umidade),
      primeiro_quartil = quantile(
        variacao_umidade,
        0.25
      ),
      terceiro_quartil = quantile(
        variacao_umidade,
        0.75
      ),
      .groups = "drop"
    ))


# Verificar o menor tamanho de grupo

resumo_resposta_precipitacao |>
  arrange(observacoes) |>
  select(
    localidade,
    variavel,
    momento,
    estado_umidade,
    classe_precipitacao,
    observacoes
  ) |>
  slice_head(n = 20)


# Localidades representativas para o texto principal

localidades_selecionadas <- c(
  "Belém",
  "Brasília",
  "Petrolina"
)


# Função para construir os gráficos

construir_grafico_precipitacao <- function(
    dados_grafico,
    momento_selecionado) {
  
  dados_grafico |>
    filter(
      momento == momento_selecionado
    ) |>
    ggplot(
      aes(
        x = classe_precipitacao,
        y = variacao_media,
        color = estado_umidade,
        group = estado_umidade
      )
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "gray50"
    ) +
    geom_line(
      linewidth = 0.7
    ) +
    geom_point(
      size = 2
    ) +
    facet_grid(
      variavel ~ localidade,
      scales = "free_y"
    ) +
    scale_x_discrete(
      labels = c(
        "Até 1 mm" = "Até 1 mm",
        "1 a 5 mm" = ">1 a 5 mm",
        "Acima de 5 mm" = "Acima de 5 mm"
      )
    ) +
    scale_color_manual(
      values = c(
        "Baixa" = "#D55E00",
        "Intermediária" = "#E69F00",
        "Alta" = "#0072B2"
      )
    ) +
    labs(
      x = "Precipitação diária",
      y = "Variação média da umidade",
      color = "Umidade anterior"
    ) +
    theme_minimal(
      base_size = 11
    ) +
    theme(
      axis.text.x = element_text(
        angle = 25,
        hjust = 1
      ),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}


# Gráficos para o texto principal

grafico_mesmo_dia_selecionadas <-
  construir_grafico_precipitacao(
    resumo_resposta_precipitacao |>
      filter(
        localidade %in% localidades_selecionadas
      ),
    "Mesmo dia"
  )

grafico_dia_seguinte_selecionadas <-
  construir_grafico_precipitacao(
    resumo_resposta_precipitacao |>
      filter(
        localidade %in% localidades_selecionadas
      ),
    "Dia seguinte"
  )


# Gráficos completos para o apêndice

grafico_mesmo_dia_completo <-
  construir_grafico_precipitacao(
    resumo_resposta_precipitacao,
    "Mesmo dia"
  )

grafico_dia_seguinte_completo <-
  construir_grafico_precipitacao(
    resumo_resposta_precipitacao,
    "Dia seguinte"
  )


# Criar as pastas

dir.create(
  "tabelas",
  showWarnings = FALSE
)

dir.create(
  "figuras",
  showWarnings = FALSE
)


# Salvar as tabelas

write_csv(
  limites_estado,
  "tabelas/limites_estados_umidade.csv"
)

write_csv(
  resumo_resposta_precipitacao,
  "tabelas/resumo_resposta_precipitacao.csv"
)


# Salvar os gráficos principais

ggsave(
  "figuras/resposta_precipitacao_mesmo_dia_selecionadas.pdf",
  grafico_mesmo_dia_selecionadas,
  width = 11,
  height = 6.5
)

ggsave(
  "figuras/resposta_precipitacao_dia_seguinte_selecionadas.pdf",
  grafico_dia_seguinte_selecionadas,
  width = 11,
  height = 6.5
)


# Salvar os gráficos completos

ggsave(
  "figuras/resposta_precipitacao_mesmo_dia_completa.pdf",
  grafico_mesmo_dia_completo,
  width = 16,
  height = 6.5
)

ggsave(
  "figuras/resposta_precipitacao_dia_seguinte_completa.pdf",
  grafico_dia_seguinte_completo,
  width = 16,
  height = 6.5
)



# A resposta da umidade à precipitação não foi uniforme entre
# as camadas, os estados anteriores e as localidades

# Para GWETTOP, chuvas acima de 5 mm produziram aumentos maiores
# quando a umidade anterior estava baixa

# Essa diferença foi mais acentuada nas localidades mais secas

# Quando a umidade anterior estava alta, a resposta foi menor,
# o que também pode estar relacionado ao limite superior da variável

# Em dias com até 1 mm de precipitação, a variação média de
# GWETTOP e GWETROOT foi negativa em todas as localidades

# GWETROOT apresentou mudanças menores que GWETTOP, mas manteve
# resposta positiva no dia seguinte às precipitações fortes

# A associação no dia seguinte possui interesse preditivo, mas
# não deve ser interpretada diretamente como efeito retardado,
# pois a precipitação do dia seguinte não foi controlada

# Os resultados indicam dependência do estado, mas não comprovam,
# isoladamente, a presença de não linearidade autorregressiva