library(INLA)
library(sf)

setwd("/dados/extremerain")
rain1=read.csv("extreme_rain_points100.csv")

data1=as.matrix(rain1$date)
idata1=order(data1)
firesco=rain1[idata1,]

library(rnaturalearth)

# -----------------------------
# ?? Brasil (robusto, sem geobr)
# -----------------------------
br <- ne_countries(country = "Brazil", returnclass = "sf")
br <- st_transform(br, 4326)
bordersp=(st_geometry(br))

fires_sf <- st_as_sf(firesco, coords = c("lon", "lat"), crs = 4326)

# filtrar pontos dentro do Brasil
fires_inside <- fires_sf[st_within(fires_sf, br, sparse = FALSE), ]

# se quiser voltar para data.frame
firesco <- cbind(st_coordinates(fires_inside), st_drop_geometry(fires_inside))


meses=as.matrix(substr(firesco$date,1,7))
umeses=unique(meses)

time=double(dim(meses)[1])+1
time[1]=1
t=1
for (j in 2:dim(meses)[1]){
  if (meses[j]!=meses[j-1]) {
    t=t+1
    time[j]=t}
  if (meses[j]==meses[j-1]) {
    time[j]=t}
}

timef=as.double(time)
obs=aggregate(rep(1,length(timef)),by=list(timef),FUN=sum)[,2]


br <- st_simplify(br, dTolerance = 0.1)

# converter para segmentos do fmesher
boundary <- fm_as_segm(br)

# pontos
locs <- firesco[, 1:2]

# mesh
smesh <- fm_mesh_2d(
  loc = locs,
  boundary = boundary,
  max.edge = c(3, 10),
  cutoff = .8
)
smesh$n
plot(smesh)



ndays <- max(time)
k <- ndays

tmesh <- fm_mesh_1d(
  loc = seq(0, ndays, length.out = k)
)


library(deldir)
dd <- deldir(smesh$loc[,1], smesh$loc[,2])
tiles <- tile.list(dd)

polys <- SpatialPolygons(lapply(1:length(tiles), function(i)
{ p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
n <- nrow(p)
Polygons(list(Polygon(p[c(1:n, 1),])), i)
}))

area <- factor(over(SpatialPoints(cbind(locs[,1], locs[,2])),
                    polys), levels=1:length(polys))

otime=time

t.breaks <- sort(c(tmesh$loc[c(1,k)],
                   tmesh$loc[2:k-1]/2 + tmesh$loc[2:k]/2))
table(time <- factor(findInterval(time, t.breaks),
                     levels=1:(length(t.breaks)-1)))

time <- factor(otime,levels=1:(length(t.breaks)-1))

agg.dat <- as.data.frame(table(area, time))
for(j in 1:2) ### set time and area as integer
  agg.dat[[j]] <- as.integer(as.character(agg.dat[[j]]))
str(agg.dat)

bounds_sf <- st_as_sf(br)

# função para calcular área de interseção
w.areas <- sapply(seq_along(tiles), function(i) {
  
  p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
  
  # fechar
  if (!all(p[1, ] == p[nrow(p), ])) {
    p <- rbind(p, p[1, ])
  }
  
  poly <- st_polygon(list(p)) |> 
    st_sfc(crs = st_crs(bounds_sf)) |>
    st_make_valid()
  
  inter <- st_intersection(poly, bounds_sf)
  
  if (length(inter) > 0) {
    as.numeric(st_area(inter))
  } else {
    0
  }
})

sum(w.areas)
summary(w.areas)

gArea_bounds <- as.numeric(st_area(bounds_sf))
gArea_bounds

w.t <- diag(fmesher::fm_fem(tmesh)$c0)

domain_area <- as.numeric(st_area(bounds_sf))


n <- nrow(firesco)

i0 <- n / (domain_area * diff(range(tmesh$loc)))
e0 <- w.areas[agg.dat$area] * (w.t[agg.dat$time])
summary(e0)

A.st <- inla.spde.make.A(
  smesh,
  loc = smesh$loc[agg.dat$area, ],
  group = agg.dat$time,
  mesh.group = tmesh
)

spde <- inla.spde2.matern(smesh)

idx <- inla.spde.make.index(
  's',
  spde$n.spde,
  n.group = k
)

loci1=smesh$loc[agg.dat$area,1:2]
vtime=agg.dat$time

rwtime=agg.dat$time
seastime=agg.dat$time
artime=agg.dat$time


stk <- inla.stack(data=list(y=agg.dat$Freq, exposure=e0),
                  A=list(A.st, 1),
                  effects=list(idx,
                               list(b0=rep(1,
                                           nrow(agg.dat)),
                                    rwtime=rwtime,
                                    seastime=seastime,
                                    artime=artime)))

save.image("modeloextremerain200.Rdata")

 

