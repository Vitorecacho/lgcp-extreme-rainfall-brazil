library(ncdf4)
library(terra)
library(sf)
library(dplyr)
library(rnaturalearth)

# -----------------------------
# ?? Brasil
# -----------------------------
br <- ne_countries(country = "Brazil", returnclass = "sf")
br <- st_transform(br, 4326)
br_vect <- vect(br)

setwd("/dados/extremerain")

# -----------------------------
# ?? Função corrigida
# -----------------------------
extract_extreme_pixels <- function(nc_file, threshold = 200) {
  
  cat("Processing:", nc_file, "\n")
  
  nc <- nc_open(nc_file)
  
  lon  <- ncvar_get(nc, "longitude")
  lat  <- ncvar_get(nc, "latitude")
  rain <- ncvar_get(nc, "precip")
  time <- ncvar_get(nc, "time")
  
  dates <- as.Date(time, origin = "1980-01-01")
  
  results <- list()
  
  for (i in seq_along(time)) {
    
    mat <- rain[,,i]
    
    # pular rápido
    if (all(mat < threshold, na.rm = TRUE)) next
    
    # criar raster
    r <- rast(
      nrows = length(lat),
      ncols = length(lon),
      xmin = min(lon), xmax = max(lon),
      ymin = min(lat), ymax = max(lat),
      crs = "EPSG:4326"
    ) 
    
    values(r) <- as.vector(t(mat))
    r <- flip(r, "vertical")
    
    # recorte Brasil
    r_br <- mask(crop(r, br_vect), br_vect)
    
    # ?? EXTRAÇÃO DIRETA (sem as.points!)
    df <- as.data.frame(r_br, xy = TRUE, na.rm = TRUE)
    
    if (nrow(df) == 0) next
    
    colnames(df)[3] <- "rain"
    
    # ?? filtro correto
    df <- df[df$rain >= threshold, ]
    
    if (nrow(df) == 0) next
    
    df$date <- dates[i]
    
    cat("  ? extremos encontrados:", nrow(df), "\n")
    
    results[[length(results) + 1]] <- df
  }
  
  nc_close(nc)
  
  if (length(results) == 0) return(NULL)
  
  bind_rows(results)
}

# -----------------------------
# ?? Rodar
# -----------------------------
files <- list.files(pattern = "\\.nc$", full.names = TRUE)

df_extremes <- lapply(files, extract_extreme_pixels) |>
  bind_rows()
  

# salvar
write.csv(df_extremes,
          "/dados/extremerain/extreme_rain_points200.csv",
          row.names = FALSE)
rm(list=ls())          
          
          