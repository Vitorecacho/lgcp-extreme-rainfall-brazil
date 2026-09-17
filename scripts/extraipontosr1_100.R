library(ncdf4)
library(terra)
library(sf)
library(dplyr)
library(rnaturalearth)

# -----------------------------
# ?? Brasil (robusto, sem geobr)
# -----------------------------
br <- ne_countries(country = "Brazil", returnclass = "sf")
br <- st_transform(br, 4326)
br_vect <- vect(br)

# -----------------------------
# ?? Diretório dos dados
# -----------------------------
setwd("/dados/extremerain")

# -----------------------------
# ?? Função principal
# -----------------------------
extract_extreme_pixels <- function(nc_file, threshold = 100) {
  
  cat("Processing:", nc_file, "\n")
  
  nc <- try(nc_open(nc_file), silent = TRUE)
  if (inherits(nc, "try-error")) return(NULL)
  
  # Variáveis CHIRPS
  lon  <- ncvar_get(nc, "longitude")
  lat  <- ncvar_get(nc, "latitude")
  rain <- ncvar_get(nc, "precip")
  time <- ncvar_get(nc, "time")
  
  # Datas (CHIRPS começa em 1980-01-01)
  dates <- as.Date(time, origin = "1980-01-01")
  
  results <- list()
  
  for (i in seq_along(time)) {
    
    mat <- rain[,,i]
    
    # pular rápido se não há extremos
    if (all(mat < threshold, na.rm = TRUE)) next
    
    # criar raster corretamente
    r <- rast(
      nrows = length(lat),
      ncols = length(lon),
      xmin = min(lon), xmax = max(lon),
      ymin = min(lat), ymax = max(lat),
      crs = "EPSG:4326"
    )
    
    values(r) <- as.vector(t(mat))
    
    # corrigir orientação (muito importante!)
    r <- flip(r, "vertical")
    
    # recorte Brasil
    r_br <- crop(r, br_vect)
    r_br <- mask(r_br, br_vect)
    
    # extrair pontos
    pts <- try(as.points(r_br, values = TRUE, na.rm = TRUE), silent = TRUE)
    if (inherits(pts, "try-error")) next
    
    vals <- values(pts)
    if (is.null(vals)) next
    
    # extremos
    pts_ext <- pts[vals >= threshold]
    if (length(pts_ext) == 0) next
    
    df <- as.data.frame(pts_ext, xy = TRUE)
    
    # segurança estrutural
    if (ncol(df) < 3) next
    
    colnames(df)[3] <- "rain"
    
    df$date <- dates[i]
    
    results[[length(results) + 1]] <- df
  }
  
  nc_close(nc)
  
  if (length(results) == 0) return(NULL)
  
  bind_rows(results)
}

# -----------------------------
# ?? Lista de arquivos
# -----------------------------
files <- list.files(pattern = "\\.nc$", full.names = TRUE)

# -----------------------------
# ? Execução (sequencial)
# -----------------------------
df_extremes <- lapply(files, extract_extreme_pixels) |>
  bind_rows()

# -----------------------------
# ?? Salvar resultado
# -----------------------------
write.csv(df_extremes, "/dados/extremerain/extreme_rain_points100.csv", row.names = FALSE)
