library(oro.nifti)
library(statcomp)
library(stringr)
library(purrr)
library(TSEntropies)
library(remotes) # Paquete usado para llamar a los repositorios de manera remota
library(neurobase)
library(wavelets)
library(planoCH)

# Función encargada de escribir una imagen nifti
writeImg <- function(img, name) {
  fname = file.path(name)
  print("Escribiendo archivo")
  writeNIfTI(img, fname, verbose = TRUE)
}

# Función para crear una imagen NIfTI vacía basada en una imagen original
create_empty_nifti <- function(original_image) {
  empty_image <- original_image
  empty_image@.Data <- array(0, dim(original_image@.Data))
  datatype(empty_image) <- 64  # FLOAT64
  bitpix(empty_image) <- 64    # 64 bits
  return(empty_image)
}

# Función para verificar si una carpeta existe y crearla si no
check_and_create_folder <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    message(paste("Carpeta creada:", folder_path))
  } else {
    message(paste("Carpeta ya existe:", folder_path))
  }
}

calculateEntropyComplexity <- function(nifti_masked_imagen, result_name, embebin) {
  # Definir los nombres de salida
  new_name_shannon_s <- str_c(result_name, "-shannon-S-ndemb", embebin)
  new_name_jsd <- str_c(result_name, "-jsd-ndemb", embebin)
  
  # Verificar si los archivos ya existen
  if (!file.exists(str_c(new_name_shannon_s, ".nii.gz")) || 
      !file.exists(str_c(new_name_jsd, ".nii.gz"))) {
    
    # Crear imágenes NIfTI vacías para los resultados
    she_image_s <- create_empty_nifti(nifti_masked_imagen)
    jsd_image <- create_empty_nifti(nifti_masked_imagen)
    
    print("Inicio de cálculos de entropía y complejidad")
    
    # Obtener dimensiones
    dims <- dim(nifti_masked_imagen)
    print(paste("Dimensiones del objeto NIfTI:", paste(dims, collapse = " x ")))
    
    if (length(dims) < 4) {
      stop("El objeto NIfTI no tiene las dimensiones esperadas (4D)")
    }
    
    # Distribuciones de referencia
    opd.unif <- rep(1 / factorial(embebin), factorial(embebin))
    opd.max <- c(1, rep(0, factorial(embebin) - 1))
    
    start_time <- Sys.time()
    
    for (z in 1:dims[3]) {  # Profundidad (z)
      start_time_z <- Sys.time()
      for (y in 1:dims[2]) {  # Altura (y)
        for (x in 1:dims[1]) {  # Anchura (x)
          if (sum(nifti_masked_imagen[x, y, z, ]) != 0) {  # Evitar voxeles vacíos
            opd <- ordinal_pattern_distribution(x = nifti_masked_imagen[x, y, z, ], ndemb = embebin)
            h <- get.shannon.disorder(opd)
            q <- get.JS.disequilibrium(opd, opd.unif, opd.max)
            
            # Asignar valores directamente a las imágenes NIfTI
            she_image_s@.Data[x, y, z, 1] <- h$S
            jsd_image@.Data[x, y, z, 1] <- q$D
          }
        }
      }
      end_time_z <- Sys.time()
      elapsed_time_z <- as.numeric(difftime(end_time_z, start_time_z, units = "secs"))
      print(paste("Tiempo para z =", z, ":", elapsed_time_z, "segundos"))
    }
    
    end_time <- Sys.time()
    total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Guardar valores constantes en los metadatos
    she_image_s@aux_file <- str_c("Smax: ", h$Smax)
    jsd_image@aux_file <- str_c("Q0: ", q$Q0)
    
    print("Fin de cálculos de entropía y complejidad")
    print("Escribiendo archivos")
    
    # Guardar resultados
    writeImg(she_image_s, new_name_shannon_s)
    writeImg(jsd_image, new_name_jsd)
    print("Archivos escritos")
    
  } else {
    print("Ya fue calculada la entropía y complejidad")
  }
}

# Ruta por defecto para la carpeta de salida
output_folder <- "result"  # Carpeta donde se guardarán los resultados

# Verificar y crear la carpeta de salida
check_and_create_folder(output_folder)

# Procesar solo los archivos proporcionados como argumentos
args <- commandArgs(trailingOnly = TRUE)

# **Mostrar los argumentos de entrada**
print("Argumentos proporcionados:")
print(args)

if (length(args) < 2) {
  stop("Proporcione al menos un archivo y el valor de 'ndemb' como argumentos.")
}

# Extraer el valor de `ndemb` del último argumento
ndemb <- as.numeric(tail(args, n = 1))
files <- head(args, n = -1)

# Validar que `ndemb` sea un número válido
if (is.na(ndemb) || ndemb <= 0) {
  stop("El valor de 'ndemb' debe ser un número positivo.")
}

# Función para procesar un archivo
process_file <- function(file_to_load, ndemb) {
  # Cargar la imagen NIfTI enmascarada
  nifti_masked_imagen <- readNIfTI(file_to_load)
  print(nifti_masked_imagen)
  
  # Extraer el nombre base del archivo para usarlo como parte del nombre de salida
  file_name <- tools::file_path_sans_ext(basename(file_to_load))
  result_name <- paste(output_folder, file_name, sep = "/")
  
  # Llama a la función para calcular entropía y complejidad
  calculateEntropyComplexity(nifti_masked_imagen, result_name, ndemb)
}

# Iterar sobre cada archivo especificado
for (file in files) {
  print(paste("Procesando archivo:", file, "con ndemb =", ndemb))
  process_file(file, ndemb)
}
