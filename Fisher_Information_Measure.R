Fisher_Information_Measure <- function(x) {
  # x: Vector numérico que representa la serie de tiempo o la distribución
  
  N <- length(x)  # Número de elementos en la serie
  if (N < 2) {
    stop("La serie debe contener al menos dos elementos.")
  }
  
  # Constante F_0, puede ajustarse según el contexto
  F_0 <- 1  # Ajustable si se requiere otra constante
  
  # Cálculo de FIM según la ecuación (35)
  fim <- F_0 * sum(((sqrt(x[1:(N-1)] + 1) - sqrt(x[1:(N-1)]))^2))
  
  return(fim)
}

# Serie de ejemplo
serie_tiempo <- c(0.1, 0.2, 0.15, 0.3, 0.25, 0.4)

# Calcular la Medida de Información de Fisher
fim_result <- Fisher_Information_Measure(serie_tiempo)
print(fim_result)
