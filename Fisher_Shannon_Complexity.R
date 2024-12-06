#Fisher_Shannon_Complexity 
# Función para calcular la Medida de Información de Fisher (FIM)
Fisher_Information_Measure <- function(x) {
  # x: Vector numérico que representa la serie de tiempo o la distribución
  
  N <- length(x)  # Número de elementos en la serie
  if (N < 2) {
    stop("La serie debe contener al menos dos elementos.")
  }
  
  # Constante F_0, puede ajustarse según el contexto
  F_0 <- 1  # Ajustable si se requiere otra constante
  
  # Cálculo de FIM según la fórmula proporcionada
  fim <- F_0 * sum(((sqrt(x[1:(N-1)] + 1) - sqrt(x[1:(N-1)]))^2))
  
  return(fim)
}

# Función para calcular el potencial de entropía de Shannon (SEP), incluyendo H(X)
SEP <- function(x) {
  # x: Vector numérico que representa la serie o la distribución
  
  # Calcular las probabilidades
  probabilidades <- table(x) / length(x)
  
  # Calcular la entropía de Shannon (H)
  H <- -sum(probabilidades * log(probabilidades))
  
  # Calcular SEP directamente
  sep <- (1 / (2 * pi * exp(1))) * exp(2 * H)
  
  return(sep)
}

# Función para calcular la Complejidad Fisher-Shannon (CFS)
Fisher_Shannon_Complexity <- function(x) {
  # x: Serie de tiempo o distribución
  
  # Calcular FIM
  fim <- Fisher_Information_Measure(x)
  
  # Calcular SEP directamente
  sep <- SEP(x)
  
  # Calcular CFS
  cfs <- fim * sep
  
  return(cfs)
}

# Ejemplo de uso
# Serie de ejemplo
serie_tiempo <- c(0.1, 0.2, 0.15, 0.3, 0.25, 0.4)

# Calcular FIM
fim_result <- Fisher_Information_Measure(serie_tiempo)
cat("Medida de Información de Fisher (FIM):", fim_result, "\n")

# Calcular SEP
sep_result <- SEP(serie_tiempo)
cat("Potencial de Entropía de Shannon (SEP):", sep_result, "\n")

# Calcular CFS (Complejidad Fisher-Shannon)
cfs_result <- Fisher_Shannon_Complexity(serie_tiempo)
cat("Complejidad Fisher-Shannon (CFS):", cfs_result, "\n")
