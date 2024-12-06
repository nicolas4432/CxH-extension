library(httr)
library(oro.nifti)

# Función para obtener el token desde variables de entorno
get_token <- function() {
  token <- Sys.getenv("GITHUB_TOKEN")
  if (token == "") stop("No se encontró el token en las variables de entorno.")
  return(token)
}

# Función para listar el contenido de una URL
listar_contenido <- function(url, token) {
  tryCatch({
    req <- GET(url, add_headers(Authorization = paste("token", token)))
    if (status_code(req) != 200) stop("Error al acceder a la URL: ", url)
    return(content(req))
  }, error = function(e) {
    cat("Error: ", e$message, "\n")
    return(NULL)
  })
}

# Función para descargar un archivo
descargar_archivo <- function(file_url, destino, token) {
  tryCatch({
    req <- GET(file_url, add_headers(Authorization = paste("token", token)))
    if (status_code(req) == 200) {
      file_content <- content(req, "raw")
      writeBin(file_content, destino)
      cat("Archivo guardado en:", destino, "\n")
    } else {
      cat("Error al descargar el archivo:", file_url, "\n")
    }
  }, error = function(e) {
    cat("Error descargando archivo:", file_url, "-", e$message, "\n")
  })
}

# Variables iniciales
token <- "inserte token aqui"

if (token == "") {
  stop("Por favor, define tu token en la variable de entorno 'GITHUB_TOKEN'.")
}
base_url <- "https://api.github.com/repos/jljara/Charite01/contents"
carpeta_descargas <- "descargas_nifti"
if (!dir.exists(carpeta_descargas)) dir.create(carpeta_descargas)

# Configuración de filtros (soporte para lista de filtros)
modo_interactivo <- FALSE
if (modo_interactivo) {
  filtro_carpeta <- readline("Ingresa el prefijo o sufijo para filtrar las carpetas: ")
  filtro_archivo <- unlist(strsplit(readline("Ingresa los formatos de archivo para filtrar, separados por comas: "), ","))
  filtro_archivo <- trimws(filtro_archivo)
  subcarpetas <- unlist(strsplit(readline("Ingresa las subcarpetas, separadas por comas: "), ","))
  subcarpetas <- trimws(subcarpetas)
} else {
  filtro_carpeta <- "sub-"
  filtro_archivo <- c("cAsym_desc-preproc_bold.nii.gz", "cAsym_desc-brain_mask.nii.gz")  # Lista de filtros
  subcarpetas <- c("func")
}

# Listar carpetas
contenido_base <- listar_contenido(base_url, token)
folders <- Filter(function(item) item$type == "dir" && grepl(filtro_carpeta, item$name), contenido_base)

if (length(folders) > 0) {
  all_file_data <- data.frame(Carpeta = character(), Nombre = character(), Tamaño = numeric(), stringsAsFactors = FALSE)
  for (folder in folders) {
    for (subcarpeta in subcarpetas) {
      subfolder_url <- paste0(base_url, "/", folder$name, "/", subcarpeta)
      contenido_subcarpeta <- listar_contenido(subfolder_url, token)
      if (!is.null(contenido_subcarpeta)) {
        # Filtra archivos que coincidan con cualquiera de los filtros
        archivos <- Filter(function(item) item$type == "file" && any(sapply(filtro_archivo, function(f) grepl(f, item$name))), contenido_subcarpeta)
        for (archivo in archivos) {
          ruta_descarga <- file.path(carpeta_descargas, archivo$name)
          file_url <- paste0("https://raw.githubusercontent.com/jljara/Charite01/main/", folder$name, "/", subcarpeta, "/", archivo$name)
          descargar_archivo(file_url, ruta_descarga, token)
          all_file_data <- rbind(all_file_data, data.frame(Carpeta = paste(folder$name, subcarpeta, sep = "/"),
                                                           Nombre = archivo$name,
                                                           Tamaño = archivo$size))
        }
      }
    }
  }
  if (nrow(all_file_data) > 0) {
    cat("\nArchivos descargados:\n")
    print(all_file_data)
  } else {
    cat("No se encontraron archivos que coincidan con el formato especificado.\n")
  }
} else {
  cat("No se encontraron carpetas que coincidan con el filtro.\n")
}
