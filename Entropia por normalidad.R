Entropia_por_normalidad <- function(serie_tiempo, alfa = 0.05, clases = 2, m = 2, base = "e", transformar = NULL, 
                                    normalizar = TRUE, probabilidades = FALSE) {
  if (!probabilidades) {
    dist_result <- Distribucion_normalidad(serie_tiempo, alfa = alfa, clases = clases, m = m, transformar = transformar)
    simbolos <- dist_result$simbolos
    probabilidades <- dist_result$probabilidades
  } else {
    simbolos <- NULL
    probabilidades <- serie_tiempo
    probabilidades <- probabilidades[probabilidades > 0]
  }
  
  funcion_log <- ifelse(base == "e", log, log2)
  
  entropia <- -sum(probabilidades * funcion_log(probabilidades))
  
  if (normalizar) {
    smax <- funcion_log(clases^m)
    entropia <- entropia / smax
  }
  
  return(list(entropia = entropia, simbolos = simbolos, probabilidades = probabilidades))
}


Distribucion_normalidad <- function(serie_tiempo, alfa = 0.05, clases = 2, m = 2, transformar = NULL) {
  patrones_normalidad <- Simbolizacion_normalidad(serie_tiempo, alfa = alfa, clases = clases, m = m, transformar = transformar)
  
  simbolos <- unique(patrones_normalidad)
  cuenta_simbolos <- table(patrones_normalidad)
  probabilidades <- as.numeric(cuenta_simbolos) / sum(cuenta_simbolos)
  
  if (round(sum(probabilidades), 10) != 1) {
    warning("Error potencial al calcular probabilidades")
  }
  
  return(list(simbolos = simbolos, probabilidades = probabilidades))
}


Simbolizacion_normalidad <- function(serie_tiempo, alfa = 0.05, clases = 2, m = 2, transformar = NULL) {
  nmuestras <- length(serie_tiempo)
  mod <- nmuestras %% m
  
  if (!is.null(transformar)) {
    serie_tiempo <- transformar(serie_tiempo)
  }
  
  sw_pval <- shapiro.test(serie_tiempo)$p.value
  
  if (sw_pval > alfa) {
    mu <- mean(serie_tiempo)
    sd <- sd(serie_tiempo)
    dist <- qnorm
    q <- seq(1 / clases, 1 - 1 / clases, length.out = clases - 1)
    b <- c(-Inf, qnorm(q, mean = mu, sd = sd), Inf)
  } else {
    q <- seq(0, 1, length.out = clases + 1)
    b <- quantile(serie_tiempo, q)
  }
  
  particiones <- cut(serie_tiempo, breaks = b, labels = FALSE, include.lowest = TRUE)
  particiones <- particiones[1:(nmuestras - mod)]
  patrones_de_normalidad <- matrix(particiones, ncol = m, byrow = TRUE)
  
  return(as.vector(patrones_de_normalidad))
}

# Serie de tiempo de ejemplo
serie_tiempo <- rnorm(100)  # Genera una serie de datos aleatorios con distribución normal

# Calcular entropía por normalidad
resultado <- Entropia_por_normalidad(serie_tiempo, alfa = 0.05, clases = 3, m = 2, base = "e")
print(resultado)
