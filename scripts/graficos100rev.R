


library(sf)
library(INLA)
setwd("/dados/extremerain")

load("/dados/extremerain/modeloextremerain100intercept.Rdata")

dates=unique(substr(firesco$date,1,7))

head(firesco)
firesco$one=1
mcount=aggregate(firesco$one,by=list(substr(firesco$date,1,7)),FUN=sum)
ts.plot(mcount[,2])
head(mcount)

library(ggplot2)
library(zoo)

# -----------------------------
# Monthly count aggregation
# -----------------------------
firesco$one <- 1

mcount <- aggregate(
  firesco$one,
  by = list(substr(firesco$date, 1, 7)),
  FUN = sum
)

colnames(mcount) <- c("month", "count")

# converter para data
mcount$time <- as.Date(
  as.yearmon(mcount$month, "%Y-%m")
)

# -----------------------------
# Plot no mesmo estilo
# -----------------------------
p_count <- ggplot(mcount, aes(x = time, y = count)) +
  
  geom_area(
    fill = "grey85",
    alpha = 0.8
  ) +
  
  geom_line(
    color = "black",
    linewidth = 0.6
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text = element_text(size = 10)
  ) +
  
  labs(
    title = "Monthly occurrence counts - 100mm"
  )

# exibir
p_count

# salvar
ggsave(
  "monthly_occurrence_counts_100.png",
  p_count,
  dpi = 600,
  width = 8,
  height = 3.5
)


par(mfrow=c(3,1))
ts.plot(resm2$summary.random$rwtime[,4:6],col=c(2,1,2))
ts.plot(resm2$summary.random$seastime[,4:6],col=c(2,1,2))
ts.plot(resm2$summary.random$artime[,4:6],col=c(2,1,2))


library(ggplot2)
library(dplyr)
library(patchwork)

dates=unique(substr(firesco$date,1,7))

# -----------------------------
# � Função para montar dataframe
# -----------------------------
make_df <- function(x, name, start_date, freq = "month") {
  
  n <- nrow(x)
  
  if (freq == "month") {
     time <- as.Date(paste0(dates, "-01"), format = "%Y-%m-%d")
  } else {
    time <- seq.Date(as.Date(start_date), by = "day", length.out = n)
  }
  
  data.frame(
    time = time,
    low = x[,4],
    median = x[,5],
    high = x[,6],
    component = name
  )
}

# -----------------------------
# � Construir datasets
# -----------------------------
df_rw   <- make_df(resm2$summary.random$rwtime,   "Trend (RW)",   "1981-01-01")
df_seas <- make_df(resm2$summary.random$seastime, "Seasonality",  "1981-01-01")
df_ar   <- make_df(resm2$summary.random$artime,   "Short-term (AR)", "1981-01-01")

# -----------------------------
# � Função de plot
# -----------------------------
plot_component <- function(df) {
  
  ggplot(df, aes(x = time)) +
    
    geom_ribbon(
      aes(ymin = low, ymax = high),
      fill = "grey80",
      alpha = 0.7
    ) +
    
    geom_line(
      aes(y = median),
      color = "black",
      linewidth = 0.6
    ) +
    
    theme_minimal(base_size = 13) +
    
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10)
    ) +
    
    labs(title = unique(df$component))
}

# -----------------------------
# � Criar plots
# -----------------------------
p1 <- plot_component(df_rw)
p2 <- plot_component(df_seas)
p3 <- plot_component(df_ar)

# -----------------------------
# � Combinar (layout Nature)
# -----------------------------
p_final <- p1 / p2 / p3

p_final

ggsave(
  "temporal_components100.png",
  p_final,
  dpi = 600,
  width = 8,
  height = 10
)


prj <- inla.mesh.projector(
  smesh,
  xlim = st_bbox(bounds_sf)[c("xmin","xmax")],
  ylim = st_bbox(bounds_sf)[c("ymin","ymax")],
  dims = c(700, 700)
)

grid_sf <- st_as_sf(
  setNames(as.data.frame(prj$lattice$loc), c("x", "y")),
  coords = c("x", "y"),
  crs = st_crs(bounds_sf)
)


inside <- st_within(grid_sf, bounds_sf, sparse = FALSE)
g.no.in <- !as.vector(inside)



n.spde <- smesh$n
k <- max(idx$s.group)
k=46

field_mat <- matrix(
  resm2$summary.random$s$mean,
  nrow = n.spde,
  ncol = k
)

t.mean <- lapply(1:k, function(j) {

  z <- inla.mesh.project(
    prj,
    field_mat[, j]   # � agora correto
  )
  
  z[g.no.in] <- NA
  
  return(z)
})



par(mfrow=c(1,4))
image(t.mean[[1]])
image(t.mean[[15]])
image(t.mean[[30]])
image(t.mean[[44]])



library(ggplot2)
library(dplyr)
library(sf)
library(viridis)
library(ggspatial)

# -----------------------------
# � Converter projeções
# -----------------------------
df_list <- lapply(seq_along(t.mean), function(j) {
  data.frame(
    x = prj$lattice$loc[,1],
    y = prj$lattice$loc[,2],
    z = as.vector(t.mean[[j]]),
    time = j,
    year = 1980 + j
  )
})

df_plot <- bind_rows(df_list) |> 
  filter(!is.na(z))

# -----------------------------
# � Selecionar anos (a cada 5)
# -----------------------------
years_sel <- seq(1981, max(df_plot$year), by = 5)

df_plot_sub <- df_plot |>
  filter(year %in% years_sel)

df_plot_sub$year <- factor(df_plot_sub$year, levels = years_sel)

# -----------------------------
# � Escala global
# -----------------------------
zlim <- range(df_plot$z, na.rm = TRUE)

# -----------------------------
# � Anotação apenas no 1º painel
# -----------------------------
annot_df <- data.frame(
  year = factor(min(years_sel), levels = levels(df_plot_sub$year))
)

# -----------------------------
# � Plot final (2x5)
# -----------------------------

zlim_all=c(-6.5,8.25)
p <- ggplot(df_plot_sub, aes(x = x, y = y, fill = z)) +
  
  geom_raster() +
  
  geom_sf(
    data = bounds_sf,
    inherit.aes = FALSE,
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  
  scale_fill_viridis_c(
    option = "magma",
    limits = zlim_all,
    name = "Spatial effect"
  ) +
  
  coord_sf(expand = FALSE) +
  
  facet_wrap(~year, nrow = 2, ncol = 5) +
  
  # � escala (apenas primeiro painel)
  annotation_scale(
    data = annot_df,
    location = "bl",
    width_hint = 0.25
  ) +
  
  # � norte (apenas primeiro painel)
 annotation_north_arrow(
  data = annot_df,
  location = "br",
  which_north = "true",
  height = unit(0.75, "cm"),
  width  = unit(0.75, "cm"),
  style = north_arrow_fancy_orienteering
) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    
    strip.text = element_text(face = "bold", size = 11),
    
    # � ajuste fino dos eixos
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10)
  ) +
  
  labs(
    title = "Estimated spatial field (SPDE)",
    subtitle = "Posterior mean of the latent spatial effect (5-year intervals)",
    x = "Longitude",
    y = "Latitude"
  )

# exibir
p

ggsave(
  "spde_field_2x5100rev.png",
  plot = p,
  dpi = 600,
  width = 12,
  height = 7,
  units = "in"
)