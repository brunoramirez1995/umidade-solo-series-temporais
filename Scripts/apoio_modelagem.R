
# Definições compartilhadas pelos scripts de modelagem

inicio_treinamento <- as.Date("2001-01-01")
fim_treinamento <- as.Date("2021-12-31")
inicio_validacao <- fim_treinamento + 1
fim_validacao <- as.Date("2024-12-31")

alvos <- c("GWETTOP", "GWETROOT")

# Para os modelos ajustados: reestimar no início de cada ano.
# As duas referências simples não têm parâmetros para reestimar.

calendario_validacao <- tibble::tibble(
  ano = 2022:2024,
  inicio_ajuste = inicio_treinamento,
  fim_ajuste = as.Date(paste0(2021:2023, "-12-31")),
  inicio_avaliacao = as.Date(paste0(2022:2024, "-01-01")),
  fim_avaliacao = as.Date(paste0(2022:2024, "-12-31"))
)

calendario_validacao


calcular_metricas <- function(observado, previsto) {
  
  stopifnot(
    length(observado) > 0,
    length(observado) == length(previsto),
    all(is.finite(observado)),
    all(is.finite(previsto))
  )
  
  erro <- observado - previsto
  
  tibble::tibble(
    n = length(erro),
    RMSE = sqrt(mean(erro^2)),
    MAE = mean(abs(erro)),
    vies = mean(erro) # Positivo: o modelo subestima, em média.
  )
}
