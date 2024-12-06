library(neurobase)

# Función para escribir una imagen NIfTI
writeImg <- function(img, name) {
  fname = file.path(name)
  cat("Escribiendo archivo en:", fname, "\n")
  writeNIfTI(img, fname, verbose = TRUE)
}

# Función para aplicar la máscara a una imagen
applyMask <- function(imagePath, maskPath) { 
  # Cargar la imagen y la máscara
  img_nifti <- readNIfTI(imagePath)
  mask_nifti <- readNIfTI(maskPath)
  
  # Extraer los datos de los objetos NIfTI
  img_data <- img_nifti@.Data
  mask_data <- mask_nifti@.Data
  
  # Verificar que las dimensiones coincidan
  if (length(dim(img_data)) == 4 && length(dim(mask_data)) == 3) {
    if (all(dim(img_data)[1:3] == dim(mask_data))) {
      # Aplica la máscara a cada volumen temporal de la imagen
      for (t in 1:dim(img_data)[4]) {
        img_data[,,,t] <- img_data[,,,t] * mask_data
      }
    } else {
      stop("Las dimensiones de la máscara no coinciden con las primeras tres dimensiones de la imagen.")
    }
  } else {
    stop("La imagen debe ser 4D (x, y, z, t) y la máscara debe ser 3D (x, y, z).")
  }
  
  # Volver a crear el objeto NIfTI con los datos modificados
  img_nifti@.Data <- img_data
  
  return(img_nifti)
}

# Función para procesar una carpeta completa
processFolder <- function(folderPath, outputFolder, suffix = "preproc_bold") {
  # Obtener una lista de archivos en el directorio
  files <- list.files(folderPath, pattern = "\\.nii\\.gz$", full.names = TRUE)
  
  # Separar máscaras e imágenes
  maskFiles <- files[grep("desc-brain_mask", files)]
  boldFiles <- files[grep(paste0("desc-", suffix), files)]
  
  # Crear la carpeta de salida si no existe
  if (!dir.exists(outputFolder)) {
    dir.create(outputFolder, recursive = TRUE)
  }
  
  # Procesar cada par de archivos
  for (boldFile in boldFiles) {
    # Identificar el sujeto y la tarea del archivo actual
    subjectTask <- gsub(paste0("_desc-", suffix, "\\.nii\\.gz$"), "", basename(boldFile))
    maskFile <- maskFiles[grep(subjectTask, maskFiles)]
    
    if (length(maskFile) == 1) {
      cat("Procesando:", boldFile, "con máscara:", maskFile, "\n")
      
      # Aplicar la máscara
      maskedImg <- tryCatch({
        applyMask(boldFile, maskFile)
      }, error = function(e) {
        cat("Error procesando:", boldFile, "\n", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(maskedImg)) {
        # Generar el nombre de salida
        outputName <- file.path(outputFolder, paste0(subjectTask, "_desc-", suffix))
        
        # Escribir la imagen procesada
        writeImg(maskedImg, outputName)
      }
    } else {
      cat("No se encontró una máscara válida para:", boldFile, "\n")
    }
  }
}

# Llamar a la función para procesar la carpeta
inputFolder <- "descargas_nifti"
outputFolder <- "procesadas_nifti"
processFolder(inputFolder, outputFolder)
