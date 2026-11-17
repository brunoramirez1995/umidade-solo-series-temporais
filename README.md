# Previsão da umidade do solo em condições climáticas contrastantes

**Autor:** Bruno Bee Ramirez  
**Disciplina:** Modelos de Séries Temporais — Universidade de São Paulo, 2026.

Este projeto compara previsões diárias de umidade do solo em seis localidades brasileiras. São avaliadas referências simples, modelos lineares, modelos com precipitação e modelos autorregressivos com regimes. Uma análise complementar examina a aplicação direta e a calibração de modelos entre localidades.

O repositório reúne os dados, scripts, tabelas, figuras e resultados necessários à reprodução das análises.

## Dados

Os dados foram extraídos do [NASA POWER](https://power.larc.nasa.gov/) em agosto de 2026, com frequência diária e tempo solar local. A base cobre 01/01/2001 a 31/12/2025 e contém 54.786 registros: 9.131 dias por localidade.

| Grupo utilizado no estudo | Localidades |
|---|---|
| Úmido | Belém e Manaus |
| Tropical sazonal | Brasília e Uberlândia |
| Mais seco | Montes Claros e Petrolina |

As variáveis-alvo são `GWETTOP` e `GWETROOT`, índices relativos de umidade na escala [0,1]. A base também contém precipitação, temperatura, umidade relativa do ar, radiação solar e evapotranspiração.

| Arquivo | Conteúdo |
|---|---|
| `Dados/nasa_power_diario.csv` | Base diária utilizada nas análises. |
| `Dados/localidades.csv` | Identificação, grupos e coordenadas das localidades. |
| `Dados/dicionario_variaveis_diarias.csv` | Descrição e unidades das variáveis. |

Para reproduzir os resultados apresentados, utilize o CSV fornecido. O script 1 documenta a coleta; uma nova consulta ao serviço pode retornar uma versão diferente dos dados.

## Protocolo de avaliação

| Etapa | Período | Uso |
|---|---|---|
| Treinamento inicial | 2001–2021 | Exploração e primeiro ajuste dos modelos. |
| Validação | 2022–2024 | Comparação e seleção das configurações. |
| Teste final | 2025 | Avaliação das configurações previamente definidas. |

As previsões são de **um dia à frente**. Na validação, os modelos ajustados são reestimados antes de cada ano, incorporando ao treinamento os dados disponíveis até dezembro do ano anterior. Os coeficientes permanecem fixos dentro do ano; o histórico observado é atualizado diariamente. As referências de persistência e sazonal ingênua não exigem estimação de parâmetros.

No teste, os ajustes utilizam dados até dezembro de 2024. As configurações são lidas de `resultados/modelagem/configuracoes_finais_2025.csv`, acompanhado por `registro_decisoes_2025.txt` na mesma pasta. Esses arquivos são entradas fornecidas: o script 20 não repete a seleção pelo desempenho em 2025.

São comparados persistência, previsão sazonal ingênua, AR, ARMA/ARIMA, AR com precipitação, AR com precipitação e interação, SETAR e LSTAR. A referência sazonal utiliza uma defasagem fixa de 365 dias; os termos harmônicos utilizam período de 365,25 dias.

Nos modelos com chuva, a previsão de Y(t+1) utiliza P(t) e, quando prevista na especificação, a interação Y(t) × P(t). A avaliação retrospectiva pressupõe a disponibilidade dessas entradas no momento da previsão.

A transferência utiliza o modelo AR com precipitação e interação em 30 direções por camada. A aplicação direta reutiliza o modelo da origem com entradas do destino. A calibração ajusta intercepto e inclinação com o histórico do destino anterior ao ano avaliado. Ambas são comparadas ao ajuste local da mesma família e à persistência.

## Organização

| Pasta | Conteúdo |
|---|---|
| `Dados/` | Base utilizada, cadastro e dicionário. |
| `Scripts/` | Scripts numerados e funções de apoio. |
| `tabelas/` | Resultados exploratórios. |
| `resultados/modelagem/` | Previsões, métricas, configurações, ajustes e diagnósticos. |
| `figuras/` | Gráficos exportados pelos scripts. |

## Ambiente e pacotes

As análises foram realizadas no **R 4.4.2, em Windows 11**. A coleta utilizou `nasapower 4.2.2`. As sessões registradas também incluem `dplyr 1.1.4`, `tidyr 1.3.1`, `readr 2.1.5` e `ggplot2 4.0.2`.

Para instalar os pacotes das análises:

```r
install.packages(c("tidyverse", "tseries"))
```

Para executar também a coleta:

```r
install.packages("nasapower")
```

Esses comandos instalam as versões disponíveis no repositório de pacotes; não fixam as versões do estudo. Os registros completos das sessões estão nos arquivos `sessionInfo.txt` das pastas de teste e transferência.

Os ajustes lineares e ARIMA utilizam funções do pacote `stats`, incluído no R. SETAR e LSTAR são implementados nos próprios scripts. Não é necessário instalar `forecast` ou `tsDyn`.

## Como executar

Execute os scripts a partir da **raiz do projeto**, onde estão `Dados/` e `Scripts/`. Todos os caminhos são relativos a essa pasta. Preserve os nomes e a capitalização das pastas utilizados pelo código.

Os arquivos `apoio_modelagem.R` e `apoio_avaliacao_final.R` são carregados pelos scripts correspondentes e devem permanecer em `Scripts/`.

Exemplo de execução no console do R:

```r
source("Scripts/2 - Verificação da Base.R", encoding = "UTF-8")
```

### Reprodução completa

| Ordem | Script(s) | Etapa |
|---|---|---|
| 1 | 2 a 11, em ordem numérica | Verificação da base, exploração, decomposição, dependência temporal, estacionariedade e covariáveis. |
| 2 | 12 | Persistência e previsão sazonal ingênua. |
| 3 | 13 | Modelos autorregressivos. |
| 4 | 14 | Modelos ARMA/ARIMA. |
| 5 | 15 | Modelos com precipitação e interação. |
| 6 | 16 | Modelos SETAR. |
| 7 | 17 | Modelos LSTAR. |
| 8 | 18 | Comparações e diagnósticos da validação. |
| 9 | 19, com `etapa <- "validacao"` | Transferência em 2022–2024. |
| 10 | 20 | Ajustes finais e avaliação local em 2025. |
| 11 | 19, com `etapa <- "teste"` | Transferência em 2025 e comparação com a validação. |

No script 19, altere a linha `etapa <- ...` **dentro do arquivo** antes de executá-lo. Definir esse objeto apenas no console não altera o modo, pois o script o redefine.

```r
# Após definir etapa <- "validacao" dentro do script 19:
source("Scripts/19 - Transferência Localidades.R", encoding = "UTF-8")

# Executar a avaliação final com as configurações fornecidas:
source("Scripts/20 - Avaliação Final.R", encoding = "UTF-8")

# Após definir etapa <- "teste" dentro do script 19:
source("Scripts/19 - Transferência Localidades.R", encoding = "UTF-8")
```

Os scripts 13–17 incluem buscas de configurações e podem exigir mais tempo. O script 20 reestima somente as configurações finais, sem repetir essas grades. Os scripts 18 e 19 reutilizam os modelos salvos; o 19 estima os coeficientes das calibrações de transferência.

A execução sobrescreve os arquivos de saída correspondentes. Execute cada script até o final para concluir as exportações e fechar os dispositivos gráficos.

### Consultar ou reutilizar resultados existentes

Os CSVs e PDFs podem ser consultados sem executar os modelos. Para reproduzir os diagnósticos do script 18 e a transferência do script 19 sem repetir os ajustes, são necessários os arquivos RDS abaixo.

Os RDS são produtos dos scripts e ficam fora do histórico do Git. Estão disponíveis no pacote [ajustes_salvos.zip](https://github.com/brunoramirez1995/umidade-solo-series-temporais/releases/download/v1.0.0/ajustes_salvos.zip), da versão v1.0.0, e também podem ser regenerados pelos scripts indicados abaixo.

| Arquivo dentro de `resultados/modelagem/` | Produzido pelo script |
|---|---|
| `autorregressivos/ajustes_selecionados.rds` | 13 |
| `arma_arima/ajustes_selecionados.rds` | 14 |
| `precipitacao/ajustes_selecionados.rds` | 15 |
| `setar/ajustes_selecionados.rds` | 16 |
| `lstar/ajustes_selecionados.rds` | 17 |
| `teste_2025/ajustes_finais.rds` | 20 |
| `teste_2025/residuos_treinamento.rds` | 20 |

Se utilizar um pacote ZIP dos ajustes, extraia-o na raiz do projeto, preservando esses caminhos. A validação do script 19 usa os ajustes de precipitação; sua etapa de teste usa os ajustes finais de 2025.

## Leitura dos resultados

- RMSE e MAE menores indicam menor erro. O viés é calculado como **observado menos previsto**: valores positivos indicam subestimação média.
- As previsões são avaliadas sem truncamento ao intervalo [0,1]. As tabelas registram as ocorrências fora dos limites.
- As colunas de ganho identificam a referência utilizada. Ganho positivo significa redução do erro em relação a essa referência.
- Resíduos de treinamento e erros de previsão são resultados distintos. Os diagnósticos dos ajustes não substituem a avaliação fora da amostra.
- Os arquivos de falhas fazem parte do registro das análises. Um CSV somente com cabeçalho indica ausência de falhas naquela etapa; tentativas de otimização malsucedidas podem aparecer em registros específicos.
- As métricas anuais utilizam o ano completo. Alguns gráficos mostram apenas um intervalo, identificado no título ou na legenda.

Os resultados locais de 2025 ficam em `resultados/modelagem/teste_2025/`. Os de transferência ficam em `resultados/modelagem/transferencia/` para validação e na subpasta `teste_2025/` para o teste.
