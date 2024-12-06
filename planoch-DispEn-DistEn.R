library(ggplot2)
library(pracma)
library(planoCH)
library(statcomp)

# Función para calcular la entropía de Shannon
# lst: Un vector de probabilidades (debe sumar 1)
# norm: Un booleano que indica si se debe normalizar la entropía
# Devuelve: El valor de la entropía (normalizada si se especifica)
get_entropy <- function(lst, norm = FALSE) {
  # Calcula la entropía de Shannon: -sum(p * log(p))
  # Se asume que lst es un vector de probabilidades
  entropy <- -sum(lst * log(lst))
  
  # Normalización opcional
  if (norm) {
    smax <- log(length(lst))
    entropy <- entropy / smax
  }
  # Retorna el valor de la entropía (normalizada o no)
  return(entropy)
}

DistEn <- function(Sig, m = 2, tau = 1, Bins = "Sturges", Logx = 2, Norm = TRUE) {
  # Validación de parámetros
  if (!is.numeric(Sig) || length(Sig) <= 10) stop("Sig debe ser un vector numérico con más de 10 elementos.")
  if (!is.integer(m) || m <= 0) stop("m debe ser un entero mayor a 0.")
  if (!is.integer(tau) || tau <= 0) stop("tau debe ser un entero mayor a 0.")
  if (!is.numeric(Logx) || Logx <= 0) stop("Logx debe ser un valor numérico positivo.")
  if (!is.logical(Norm)) stop("Norm debe ser un valor lógico.")
  if (!(is.numeric(Bins) && Bins > 1) && !(tolower(Bins) %in% c("sturges", "sqrt", "rice", "doanes"))) {
    stop("Bins debe ser un número mayor a 1 o un método válido: 'sturges', 'sqrt', 'rice', 'doanes'.")
  }
  
  Sig <- as.vector(Sig)
  Nx <- length(Sig) - (m - 1) * tau
  
  # Reconstrucción del espacio
  Zm <- matrix(0, nrow = Nx, ncol = m)

  for (n in 1:m) {
    Zm[, n] <- Sig[(1:Nx) + (n - 1) * tau]
  }

  # Cálculo de distancias máximas
  DistMat <- numeric(choose(Nx, 2))
  
  
  idx <- 1
  for (k in 1:(Nx - 1)) {
    DistMat[idx:(idx + Nx - k - 1)] <- apply(abs(sweep(Zm[(k + 1):Nx, , drop = FALSE], 2, Zm[k, ])), 1, max)
    idx <- idx + (Nx - k)
  }
  
  Ny <- length(DistMat)
  
  # Selección de bins
  if (is.character(Bins)) {
    Bx <- switch(tolower(Bins),
                 sturges = ceiling(log2(Ny) + 1),
                 sqrt = ceiling(sqrt(Ny)),
                 rice = ceiling(2 * Ny^(1/3)),
                 doanes = {
                   sigma <- sqrt(6 * (Ny - 2) / ((Ny + 1) * (Ny + 3)))
                   ceiling(1 + log2(Ny) + log2(1 + abs(mean(DistMat)/sigma)))
                 },
                 stop("Método de binning no válido."))
  } else {
    Bx <- as.integer(Bins)
  }
  
  # Ajustar límites de bins manualmente
  bin_edges <- seq(min(DistMat), max(DistMat), length.out = 18)
  
  # Crear el histograma con límites ajustados
  hist_data <- hist(DistMat, breaks = bin_edges, plot = FALSE)
  Ppi <- hist_data$density / sum(hist_data$density)
  

  
  
  
  # Cálculo de entropía
  Ppi <- Ppi[Ppi > 0]
  Dist <- -sum(Ppi * log(Ppi) / log(Logx))

  
  if (Norm) {
    Dist <- Dist / (log(Bx) / log(Logx))
  }
  
  return(list(Dist = Dist, Ppi = Ppi))
}

DispEn <- function(Sig, m = 2, tau = 1, c = 3, Typex = "NCDF", Logx = exp(1), Fluct = FALSE, Norm = FALSE, rho = 1) {
  # Validación de parámetros
  if (!is.numeric(Sig) || length(Sig) <= 10) stop("Sig debe ser un vector numérico con más de 10 elementos.")
  if (!is.integer(m) || m <= 0) stop("m debe ser un entero mayor a 0.")
  if (!is.integer(tau) || tau <= 0) stop("tau debe ser un entero mayor a 0.")
  if (!is.integer(c) || c <= 1) stop("c debe ser un entero mayor a 1.")
  if (!is.numeric(Logx) || Logx <= 0) stop("Logx debe ser un valor positivo.")
  if (!is.logical(Fluct)) stop("Fluct debe ser un valor lógico.")
  if (!is.logical(Norm)) stop("Norm debe ser un valor lógico.")
  if (!Typex %in% c("linear", "kmeans", "ncdf", "finesort", "equal")) {
    stop("Typex debe ser uno de: 'linear', 'kmeans', 'ncdf', 'finesort', 'equal'.")
  }
  
  N <- length(Sig)
  Sig <- as.vector(Sig)
  
  # Transformación según Typex
  if (tolower(Typex) == "linear") {
    breaks <- seq(min(Sig), max(Sig), length.out = c + 1)
    Zi <- findInterval(Sig, breaks, rightmost.closed = TRUE)
    # cat("Depuración Zi (Typex = linear, R):\n")
    # print(head(Zi, 10))
    # print(tail(Zi, 10))
    
    } else if (tolower(Typex) == "kmeans") {
      set.seed(42)  # Semilla para reproducibilidad
      initial_centroids <- seq(min(Sig), max(Sig), length.out = c)  # Inicialización uniforme
      kmeans_res <- kmeans(Sig, centers = initial_centroids, iter.max = 200)
      
      # Ordenar centroides y reasignar etiquetas
      ordered_centroids <- order(kmeans_res$centers)
      Zi <- as.integer(factor(kmeans_res$cluster, levels = ordered_centroids))
      
      # cat("Depuración Zi (Typex = kmeans, R):\n")
      # print(head(Zi, 10))
      # print(tail(Zi, 10))
    } else if (tolower(Typex) == "ncdf") {
    Zx <- pnorm((Sig - mean(Sig)) / sqrt(sum((Sig - mean(Sig))^2) / (length(Sig) - 1)))
    Zi <- cut(Zx, breaks = seq(0, 1, length.out = c + 1), labels = FALSE, include.lowest = TRUE, right = FALSE)
    # cat("Depuración Zi (Typex = ncdf, R):\n")
    # print(head(Zi, 10))
    # print(tail(Zi, 10))
    
    } else if (tolower(Typex) == "finesort") {
      if (rho <= 0) stop("rho debe ser mayor a 0 para 'finesort'.")
      
      Zx <- pnorm((Sig - mean(Sig)) / sd(Sig))
      Zi <- findInterval(Zx, seq(0, 1, length.out = c + 1), rightmost.closed = TRUE)
      
      Ym <- matrix(NA, nrow = N - (m - 1) * tau, ncol = m)
      for (n in 1:m) {
        Ym[, n] <- Zx[(1:(N - (m - 1) * tau)) + (n - 1) * tau]
      }
      
      diff_Ym <- abs(apply(Ym, 1, diff))
      Yi <- floor(apply(diff_Ym, 1, max) / (rho * sd(abs(diff(Sig)))))
      
      # Agregar Yi a Zm
      Zm <- cbind(Ym, Yi)
      
      # cat("Depuración Zi (Typex = finesort, R):\n")
      # print(head(Zi, 10))
      # print(tail(Zi, 10))
    }
  else if (tolower(Typex) == "equal") {
    idx <- order(Sig, decreasing = FALSE)  # Ordenar los índices
    intervals <- seq(1, N, length.out = c + 1)  # Definir los intervalos
    Zi <- rep(NA, N)  # Inicializar el vector de símbolos
    
    for (i in 1:c) {
      # Asignar símbolos a los intervalos
      Zi[idx[ceiling(intervals[i]):floor(intervals[i + 1])]] <- i
    }
    
    # Reemplazar posibles NA con el valor más cercano (opcional)
    if (any(is.na(Zi))) {
      Zi[is.na(Zi)] <- c  # Asignar el último símbolo
    }
    
    # cat("Depuración Zi (Typex = equal, R):\n")
    # print(head(Zi, 10))
    # print(tail(Zi, 10))
  }
  
  # Mensaje de depuración para Zx en R
   # cat("Depuración Zx (R):\n")
   # cat(sprintf("Tamaño: %d\n", length(Zx)))
   # cat("Primeros 10 valores:\n")
   # print(head(Zx, 10))
   # cat("Últimos 10 valores:\n")
   # print(tail(Zx, 10))
   # cat(sprintf("Rango: Min=%.6f, Max=%.6f\n", min(Zx), max(Zx)))
  
  
  # Reconstrucción del espacio embebido
  Zm <- matrix(NA, nrow = N - (m - 1) * tau, ncol = m)
  for (n in 1:m) {
    Zm[, n] <- Zi[(1:(N - (m - 1) * tau)) + (n - 1) * tau]
  }


  if (tolower(Typex) == "finesort") {
    Zm <- cbind(Zm, Yi)
  }
  
  if (Fluct) {
    Zm <- diff(Zm, differences = 1)
    if (m < 2) warning("Fluctuation-based Dispersion Entropy no está definida para m = 1.")
  }
  
  # #Mensaje de depuración para Zm en R
  # cat("Depuración Zm (R):\n")
  # cat(sprintf("Dimensiones: %d x %d\n", nrow(Zm), ncol(Zm)))
  # cat("Primeros 5 vectores:\n")
  # print(head(Zm, 5))
  # cat("Últimos 5 vectores:\n")
  # print(tail(Zm, 5))

  
  # Ordenar los vectores únicos para que coincidan con Python
  unique_rows <- unique(as.data.frame(Zm))  # Asegúrate de que Zm esté como data.frame
  unique_rows <- unique_rows[do.call(order, unique_rows), ]  # Ordenar por columnas
  
  # Calcular las frecuencias de los vectores únicos
  counts <- table(do.call(paste, as.data.frame(Zm)))  # Tabla de frecuencias
  Ppi <- as.numeric(counts) / nrow(Zm)  # Probabilidades
  
  # Depuración de vectores únicos y frecuencias
  # cat("Vectores únicos ordenados:\n")
  # print(unique_rows)
  # cat("Frecuencias (conteos):\n")
  # print(as.numeric(counts))
  # cat("Probabilidades (Ppi):\n")
  # print(Ppi)
  # cat(sprintf("Suma de probabilidades: %.6f\n", sum(Ppi)))
  
  
  
  Ppi <- Ppi[Ppi > 0]  # Elimina las probabilidades iguales a 0
  
  # cat("Depuración Ppi (R):\n")
  # print(Ppi)
  # cat(sprintf("Suma: %.6f\n", sum(Ppi)))
  
  Dispx <- -sum(Ppi * log(Ppi) / log(Logx))
  
  
  if (Norm) {
    Nx <- nrow(unique_rows)
    if (Fluct) {
      Dispx <- Dispx / (log((2 * c - 1)^(m - 1)) / log(Logx))
    } else {
      Dispx <- Dispx / (log(c^m) / log(Logx))
    }
  }
  return(list(Dispx = Dispx, Ppi = Ppi))
}


# Función para calcular la constante de normalización
# n: El tamaño del conjunto (longitud del vector de probabilidades)
# dist: Un string que indica el tipo de distancia ('euc', 'woo', 'sjd')
# Devuelve: El valor de la constante de normalización basado en el tipo de distancia
get_norm_constant <- function(n, dist = 'euc') {
  q0 <- NULL
    # Caso 1: Distancia Euclidiana ('euc')
  if (dist == 'euc') { 
    q0 <- n / (n - 1)
    # Caso 2: Distancia de Wootters ('woo')
  } else if (dist == 'woo') {
    q0 <- 1 / acos((1 / n)^(1 / 2))
    # Caso 3: Divergencia de Shannon-Jensen ('sjd')
  } else if (dist == 'sjd') {
    q0 <- 1 / log(n)
  }
  # Retorna el valor de la constante de normalización q0
  return(q0)
}


# Función para calcular la distancia euclidiana
# p1: Un vector de probabilidades
get_euclidean_distance <- function(p1) {
  return(sum((p1 - (1 / length(p1)))^2))
}


# Función para calcular la distancia de Wooters
# p1: Un vector de probabilidades
get_wooters_distance <- function(p1) {
  return(acos(sum(sqrt(p1) * sqrt(1 / length(p1)))))
}


# Función para calcular la divergencia de Shannon-Jensen (SJD)
# p1: Un vector de probabilidades
# Devuelve: La divergencia de Shannon-Jensen normalizada
get_shannon_jensen_divergence <- function(p1) {
  uniform_dist <- rep(1 / length(p1), length(p1))
  p_plus_u_over_2 <- (uniform_dist + p1) / 2
  s_of_p_plus_u_over_2 <- -sum(p_plus_u_over_2 * log(p_plus_u_over_2))
  probabilities <- p1[p1 != 0]
  s_of_p_over_2 <- -sum(probabilities * log(probabilities)) / 2
  s_of_u_over_2 <- log(length(p1)) / 2
  
  # Calcula el valor máximo de la divergencia de Shannon-Jensen
  # Esto asegura que la divergencia esté normalizada
  js_div_max <- -0.5 * (((length(p1) + 1) / length(p1)) * log(length(p1) + 1) + 
                          log(length(p1)) - 2 * log(2 * length(p1)))
  
  # Calcula la divergencia de Shannon-Jensen como la diferencia de las entropías
  js_div <- s_of_p_plus_u_over_2 - s_of_p_over_2 - s_of_u_over_2
  
  # Devuelve la divergencia de Shannon-Jensen normalizada
  return(js_div / js_div_max)
}


# Función para calcular la divergencia de Kullback-Leibler (KL)
# p: Un vector de probabilidades (distribución P)
# q: Un vector de probabilidades (distribución Q)
# Devuelve: El valor de la divergencia de Kullback-Leibler (KL) entre P y Q
kl_divergence <- function(p, q) {
  return(sum(p * log2(p / q)))
}


# Función para calcular el desequilibrio
# p: Un vector de probabilidades
# dist: El tipo de distancia a usar ('euc', 'woo', 'sjd')
# Devuelve: El valor del desequilibrio basado en el tipo de distancia
get_disequilibrium <- function(p, dist = 'euc') {
  n <- length(p)  # Longitud del vector
  q0 <- get_norm_constant(n, dist)  # Constante de normalización
  
  # Calcula la distancia en función del tipo especificado
  if (dist == 'euc') {
    prob_distance <- get_euclidean_distance(p)  # Distancia euclidiana
    diseq <- q0 * prob_distance
  } else if (dist == 'woo') {
    prob_distance <- get_wooters_distance(p)  # Distancia de Wooters
    diseq <- q0 * prob_distance
  } else if (dist == 'sjd') {
    diseq <- get_shannon_jensen_divergence(p)  # Divergencia de Shannon-Jensen
  }
  return(diseq)
}


# Función para generar una distribución equiprobable
# n: El tamaño de la distribución (número de elementos)
# Devuelve: Un vector de tamaño n donde cada elemento es 1/n (distribución uniforme)
get_ref_prob_dist <- function(n) {
  return(rep(1 / n, n))  # Genera una distribución uniforme de tamaño n
}


# Función para calcular la complejidad estadística
# p: Un vector de probabilidades
# h: El valor de la entropía
# dist: El tipo de distancia a usar ('euc', 'woo', 'sjd')
# Devuelve: El valor de la complejidad estadística
get_complexity <- function(p, h, dist) {
  q <- get_disequilibrium(p, dist)  # Calcula el desequilibrio
  return(h * q)  # Multiplica el desequilibrio por la entropía para obtener la complejidad
}


# Función para calcular la entropía de permutación
# probabilities: Vector de probabilidades
# dx: Dimensión de incrustación (embedding dimension)
# dy: Dimensión de incrustación vertical (en series temporales se usa 1)
# Devuelve: La entropía de permutación
permutation_entropy <- function(probabilities, dx, dy) {
  -sum(probabilities * log(probabilities), na.rm = TRUE)
}


# Función para calcular la curva mínima de complejidad y entropía
# dx: Dimensión de incrustación (embedding dimension) horizontal
# dy: Dimensión de incrustación vertical (1 para series temporales)
# size: Longitud del vector de datos devueltos (pares de entropía y complejidad)
# Devuelve: Matriz de valores de entropía y complejidad que delimitan la curva inferior del plano de causalidad
# Función para calcular la curva mínima de complejidad y entropía
minimum_complexity_entropy <- function(dx = 3, dy = 1, size = 100) {
  size <- size + 1  # Incrementa el tamaño para incluir el valor superior
  N <- factorial(dx * dy)  # Número de estados posibles
  prob_params <- seq(1 / N, 1, length.out = size - 1)  # Genera parámetros de probabilidad
  uniform_dist <- rep(1 / N, N)  # Distribución uniforme
  
  hc <- matrix(nrow = size - 1, ncol = 2)  # Inicializa la matriz para almacenar entropía y complejidad
  
  for (i in seq_len(size - 1)) {
    probabilities <- rep((1 - prob_params[i]) / (N - 1), N)  # Inicializa las probabilidades
    probabilities[1] <- prob_params[i]  # Ajusta la primera probabilidad
    
    # Calcula la entropía de permutación
    h <- permutation_entropy(probabilities, dx, dy) / log(N)  # Normaliza la entropía con log(N)
    
    # Calcula la divergencia de Shannon-Jensen
    js_div <- get_shannon_jensen_divergence(probabilities)
    
    hc[i, ] <- c(h, h * js_div)  # Añade el par de entropía y complejidad
  }
  
  return(hc)
}


# Función para calcular la curva de complejidad y entropía máxima
maximum_complexity_entropy <- function(dx = 3, dy = 1, m = 1) {
  N <- factorial(dx * dy)
  hlist <- numeric()
  clist <- numeric()
  
  for (i in seq_len(N - 1)) {
    p <- rep(0, N)
    uniform_dist <- rep(1 / N, N)
    prob_params <- seq(0, 1 / N, length.out = m)
    
    for (k in seq_len(m)) {
      # Ajustar la primera probabilidad y luego las restantes como en Python
      p[1] <- prob_params[k]  
      p[2:(N - i)] <- (1 - prob_params[k]) / (N - i - 1)  # Ajustar correctamente las demás posiciones
      
      # Calcular la entropía de permutación y la divergencia Shannon-Jensen
      h <- permutation_entropy(p, dx, dy) / log(N)
      js_div <- get_shannon_jensen_divergence(p)
      
      # Verificar si h o js_div son NaN o Inf, si no, asignar los valores
      if (!is.nan(h) && !is.nan(js_div) && is.finite(h) && is.finite(js_div)) {
        hlist <- c(hlist, h)
        clist <- c(clist, h * js_div)
      }
    }
  }
  
  
  args <- order(hlist)
  
  return(cbind(hlist[args], clist[args]))
}



###################### Complejidad y Entropia de la señal Stocastica ######################
dx <- 3


signal_data <- read.table("stochastic_signal.txt", sep = ",", header = FALSE)

# Extraer la columna V1 como un vector numérico
sig_formatted <- signal_data$V1

# Verificar el formato del nuevo vector
#sig_formatted

# Distribución de patrones ordinales, de dimension (ndem) = 5
opd = ordinal_pattern_distribution(x = sig_formatted, ndemb = dx)

# Calculo de la entropía de permutación para comparacion
#permutation_entropy(opd)

# Definicion de distribuciones de referencia
# Crea una distribución uniforme donde cada patrón ordinal tiene la misma 
# probabilidad. Representa el caso de máxima entropía (máximo desorden).
opd.unif <- rep(1,length(opd)) 

#Crea una distribución donde solo un patrón tiene probabilidad 1 y los demás 0. 
#Representa el caso de mínima entropía (máximo orden).
opd.max <- c(1,rep(0,length(opd)-1))

#Ecuacion 45 cap 8
#Cálculo de la Complejidad Estadística usando Entropía de Shannon y Desequilibrio de Jensen-Shannon
h <- get.shannon.disorder(opd) #Calcula la entropía de Shannon                
q <- get.JS.disequilibrium(opd, opd.unif, opd.max  ) #Calcula el desequilibrio de Jensen-Shannon 
C.s.js <- h$H*q$Q #kapa = shannon = "S", nu = jensen = "Js", q = 1
#print(paste("h ", h, "q",q))
#cat("C[k=s, v=j, q=1] = ", C.s.js, "\n")

# Calcular la entropía de Shannon de la señal
entropy_signal <- h$H
#entropy_signal2 <- get_entropy(sig_formatted)

# Calcular la complejidad estadística de la señal (usando distancia Euclidiana)
complexity_signal <- C.s.js


# Cargar las señales
oscillation_signal_data <- read.table("oscillation_signal.txt", sep = ",", header = FALSE)
periodic_signal_data <- read.table("periodic_signal.txt", sep = ",", header = FALSE)

# Formatear las señales
oscillation_formatted <- oscillation_signal_data$V1
periodic_formatted <- periodic_signal_data$V1

# Calcular los patrones ordinales
opd_oscillation <- ordinal_pattern_distribution(x = oscillation_formatted, ndemb = dx)
opd_periodic <- ordinal_pattern_distribution(x = periodic_formatted, ndemb = dx)

# Calcular la entropía y complejidad para la señal de oscilación
h_oscillation <- get.shannon.disorder(opd_oscillation)
q_oscillation <- get.JS.disequilibrium(opd_oscillation, opd.unif, opd.max)
C_oscillation <- h_oscillation$H * q_oscillation$Q

# Calcular la entropía y complejidad para la señal periódica
h_periodic <- get.shannon.disorder(opd_periodic)
q_periodic <- get.JS.disequilibrium(opd_periodic, opd.unif, opd.max)
C_periodic <- h_periodic$H * q_periodic$Q

# Crear data frames con los puntos de las nuevas señales
df_oscillation_point <- data.frame(H = h_oscillation$H, C = C_oscillation, Curva = "Señal Oscilatoria")
df_periodic_point <- data.frame(H = h_periodic$H, C = C_periodic, Curva = "Señal Periódica")

# Unir los puntos de las nuevas señales con los datos existentes
#df_combined <- rbind(df_combined, df_oscillation_point, df_periodic_point)




#### PLANO C/H
# Generar las curvas mínima y máxima de complejidad-entropía
hc_min_curve <- minimum_complexity_entropy(dx = dx, size = 100)
hc_max_curve <- maximum_complexity_entropy(dx = dx, m = 2)

# Convertir las curvas en data frames para ggplot2
df_min <- data.frame(H = hc_min_curve[,1], C = hc_min_curve[,2], Curva = "Curva Mínima")
df_max <- data.frame(H = hc_max_curve[,1], C = hc_max_curve[,2], Curva = "Curva Máxima")

# Unir ambas curvas en un solo data frame
df_combined <- rbind(df_min, df_max)

# Crear un data frame con el punto de la señal estocástica
df_signal_point <- data.frame(H = entropy_signal, C = complexity_signal, Curva = "Señal Estocástica")

# Unir el nuevo punto con los datos existentes de las curvas mínima y máxima
df_combined <- rbind(df_combined, df_signal_point, df_oscillation_point, df_periodic_point)

# Graficar las curvas con los nuevos puntos añadidos
plot2 <- ggplot(df_combined, aes(x = H, y = C, color = Curva)) +
  geom_line(size = 1.5) +  # Curvas mínima y máxima
  geom_point(data = df_signal_point, aes(x = H, y = C), color = "green", size = 4) +  # Punto de la señal estocástica
  geom_point(data = df_oscillation_point, aes(x = H, y = C), color = "orange", size = 4) +  # Punto de la señal oscilatoria
  geom_point(data = df_periodic_point, aes(x = H, y = C), color = "purple", size = 4) +  # Punto de la señal periódica
  labs(x = 'Permutation Entropy (H)', y = 'Statistical Complexity (C)', title = 'Curvas Mínima, Máxima y Señales') +
  scale_color_manual(values = c("Curva Mínima" = "blue", "Curva Máxima" = "red", 
                                "Señal Estocástica" = "green", "Señal Oscilatoria" = "orange", 
                                "Señal Periódica" = "purple")) +
  theme_minimal()
print(plot2)



get.opd <- function(N = 296) {
  opd <- rep(0, 120)
  opd[1] <- sample(70:80, 1)
  
  w <- 1 / (1.5 ** (1:12))# (2*(1:12))
  w <- w / sum(w)
  
  i <- sample(2:119, 120 - opd[1])
  opd[i] <- sample(1:12, length(i), replace = TRUE, prob = w)
  opd[120] <- N - sum(opd[-120])
  
  opd
}



# Cambiar la dimensión de embebido, el retardo y el método de binning
#resultado <- DistEn(oscillation_formatted, m = as.integer(2), tau = as.integer(1), Bins = "sturges")

# Mostrar los resultados
#print(resultado$Dist)  # Entropía sin normalizar
#print(resultado$Ppi)   # Probabilidades

# Cargar la señal desde un archivo o generar una de ejemplo
#signal <- scan("oscillation_signal.txt", sep = ",")  # Reemplaza con tu archivo

# Cargar la señal desde un archivo
signal <- scan("oscillation_signal.txt", sep = ",")

# Verificar que la señal no esté vacía y sea válida
if (!is.numeric(signal) || length(signal) <= 10) {
  stop("La señal debe ser un vector numérico con más de 10 elementos.")
}

# Lista de métodos
methods <- c("linear", "kmeans", "ncdf", "finesort", "equal")

# Parámetros comunes
m <- as.integer(3)  # Asegurar que sea un entero
tau <- as.integer(2)
c <- as.integer(4)
Logx <- exp(1)  # Cambiar si necesitas base logarítmica distinta

# Lista para almacenar los resultados
results <- list()

# Loop para calcular la entropía de dispersión para cada método
for (method in methods) {
  tryCatch({
    result <- DispEn(Sig = signal, m = m, tau = tau, c = c, Typex = method, Logx = Logx, Norm = TRUE)
    results[[method]] <- result$Dispx
    cat(sprintf("Método: %s, Entropía: %f\n", method, result$Dispx))
  }, error = function(e) {
    cat(sprintf("Error con el método %s: %s\n", method, e$message))
  })
}

# Resultados finales
cat("\nResultados finales:\n")
for (method in names(results)) {
  cat(sprintf("%s: %f\n", method, results[[method]]))
}

