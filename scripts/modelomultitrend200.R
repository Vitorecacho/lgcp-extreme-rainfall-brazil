library(INLA)
library(sf)
library(fmesher)
library(sp)
setwd("/dados/extremerain")


load("modeloextremerain200.Rdata")


dates=unique(substr(firesco$date,1,7))
udates=unique(dates)
idates=1:length(udates)
gwrtime=as.double(substr(udates[rwtime],1,4))-1980
max(gwrtime)


library(zoo)

ym <- as.yearmon(udates[rwtime], "%Y-%m")
head(ym)

# índice mensal contínuo
idx_month <- as.integer(12 * (floor(ym) - min(floor(ym))) +
                                cycle(ym))

idx_month <- idx_month - 1
                  
seastime=idx_month



A.st <- inla.spde.make.A(
  smesh,
  loc = smesh$loc[agg.dat$area, ],
  group = gwrtime,
  mesh.group = tmesh
)

spde <- inla.spde2.matern(smesh)

idx <- inla.spde.make.index(
  's',
  spde$n.spde,
  n.group = max(gwrtime)
)

library(sf)
library(dplyr)
library(geobr)

# ============================================
# Estados Brasil
# ============================================

br_states <- read_state(year = 2020)

# verificar CRS
st_crs(br_states)

# ============================================
# Nós do mesh
# ============================================

mesh_coords <- data.frame(
  lon = smesh$loc[,1],
  lat = smesh$loc[,2]
)

# ============================================
# Converter para sf
# ============================================

mesh_sf <- st_as_sf(
  mesh_coords,
  coords = c("lon","lat"),
  crs = 4326
)

# ============================================
# Ajustar CRS
# ============================================

mesh_sf <- st_transform(mesh_sf, st_crs(br_states))

# ============================================
# Spatial join
# ============================================

mesh_join <- st_join(
  mesh_sf,
  br_states[,c("abbrev_state","name_region")]
)

# ============================================
# Região dos nós
# ============================================

mesh_regions <- mesh_join$name_region

table(mesh_regions)

mesh_regions[is.na(mesh_regions)]="Fora"
Si=rep(mesh_regions=="Sul",518)
Ni=rep(mesh_regions=="Norte",518)
NEi=rep(mesh_regions=="Nordeste",518)
SEi=rep(mesh_regions=="Sudeste",518)
COi=rep(mesh_regions=="Centro Oeste",518)

Srwtime=rwtime
Nrwtime=rwtime
NErwtime=rwtime
SErwtime=rwtime
COrwtime=rwtime

Srwtime=rwtime
Nrwtime=rwtime
NErwtime=rwtime
SErwtime=rwtime
COrwtime=rwtime

Srwtime=rwtime
Nrwtime=rwtime
NErwtime=rwtime
SErwtime=rwtime
COrwtime=rwtime

Sseastime=seastime
Nseastime=seastime
NEseastime=seastime
SEseastime=seastime
COseastime=seastime

Sartime=artime
Nartime=artime
NEartime=artime
SEartime=artime
COartime=artime

Srwtime[!Si]=NA
Nrwtime[!Ni]=NA
NErwtime[!NEi]=NA
SErwtime[!SEi]=NA
COrwtime[!COi]=NA

Sseastime[!Si]=NA
Nseastime[!Ni]=NA
NEseastime[!NEi]=NA
SEseastime[!SEi]=NA
COseastime[!COi]=NA

Sartime[!Si]=NA
Nartime[!Ni]=NA
NEartime[!NEi]=NA
SEartime[!SEi]=NA
COartime[!COi]=NA

stk <- inla.stack(
  data = list(
    y = agg.dat$Freq,
    exposure = e0
  ),
  
  A = list(A.st, 1),
  
  effects = list(
    
    idx,
    
    list(
      Intercept = rep(1, nrow(agg.dat)),
      
      rwtime   = seastime,
      seastime = seastime,
      artime   = seastime,
      
      Srwtime  = Srwtime,
      Nrwtime  = Nrwtime,
      NErwtime = NErwtime,
      SErwtime = SErwtime,
      COrwtime = COrwtime,
      
      Sseastime  = Sseastime,
      Nseastime  = Nseastime,
      NEseastime = NEseastime,
      SEseastime = SEseastime,
      COseastime = COseastime,
      
      Sartime  = Sartime,
      Nartime  = Nartime,
      NEartime = NEartime,
      SEartime = SEartime,
      COartime = COartime
    )
  )
)


model2r <- y ~ -1 + Intercept +

  # =====================================================
  # RW2 por região
  # =====================================================

  f(Srwtime,
    model = "rw2",
    constr = TRUE,
    scale.model = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(10, 0.1),
        initial = log(118.253),
        fixed = FALSE
      )
    )
  ) +

  f(Nrwtime,
    model = "rw2",
    constr = TRUE,
    scale.model = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(10, 0.1),
        initial = log(118.253),
        fixed = FALSE
      )
    )
  ) +

  f(NErwtime,
    model = "rw2",
    constr = TRUE,
    scale.model = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(10, 0.1),
        initial = log(118.253),
        fixed = FALSE
      )
    )
  ) +

  f(SErwtime,
    model = "rw2",
    constr = TRUE,
    scale.model = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(10, 0.1),
        initial = log(118.253),
        fixed = FALSE
      )
    )
  ) +

  f(COrwtime,
    model = "rw2",
    constr = TRUE,
    scale.model = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(10, 0.1),
        initial = log(118.253),
        fixed = FALSE
      )
    )
  ) +

  # =====================================================
  # Sazonalidade por região
  # =====================================================

  f(Sseastime,
    model = "seasonal",
    season.length = 12,
    constr = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(1, 0.01),
        initial = log(278.591),
        fixed = FALSE
      )
    )
  ) +

  f(Nseastime,
    model = "seasonal",
    season.length = 12,
    constr = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(1, 0.01),
        initial = log(278.591),
        fixed = FALSE
      )
    )
  ) +

  f(NEseastime,
    model = "seasonal",
    season.length = 12,
    constr = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(1, 0.01),
        initial = log(278.591),
        fixed = FALSE
      )
    )
  ) +

  f(SEseastime,
    model = "seasonal",
    season.length = 12,
    constr = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(1, 0.01),
        initial = log(278.591),
        fixed = FALSE
      )
    )
  ) +

  f(COseastime,
    model = "seasonal",
    season.length = 12,
    constr = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(1, 0.01),
        initial = log(278.591),
        fixed = FALSE
      )
    )
  ) +

  # =====================================================
  # AR(2) por região
  # =====================================================

  f(Sartime,
    model = "ar",
    order = 2,
    constr = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(1, 0.01),
        initial = log(1.118),
        fixed = FALSE
      ),
      pacf1 = list(
        prior = "normal",
        param = c(0, 0.15),
        initial = atanh(0.047),
        fixed = FALSE
      ),
      pacf2 = list(
        prior = "normal",
        param = c(0, 0.15),
        initial = atanh(0.056),
        fixed = FALSE
      )
    )
  ) +

  f(Nartime,
    model = "ar",
    order = 2,
    constr = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(1, 0.01),
        initial = log(1.118),
        fixed = FALSE
      ),
      pacf1 = list(
        prior = "normal",
        param = c(0, 0.15),
        initial = atanh(0.047),
        fixed = FALSE
      ),
      pacf2 = list(
        prior = "normal",
        param = c(0, 0.15),
        initial = atanh(0.056),
        fixed = FALSE
      )
    )
  ) +

  f(NEartime,
    model = "ar",
    order = 2,
    constr = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(1, 0.01),
        initial = log(1.118),
        fixed = FALSE
      ),
      pacf1 = list(
        prior = "normal",
        param = c(0, 0.15),
        initial = atanh(0.047),
        fixed = FALSE
      ),
      pacf2 = list(
        prior = "normal",
        param = c(0, 0.15),
        initial = atanh(0.056),
        fixed = FALSE
      )
    )
  ) +

  f(SEartime,
    model = "ar",
    order = 2,
    constr = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(1, 0.01),
        initial = log(1.118),
        fixed = FALSE
      ),
      pacf1 = list(
        prior = "normal",
        param = c(0, 0.15),
        initial = atanh(0.047),
        fixed = FALSE
      ),
      pacf2 = list(
        prior = "normal",
        param = c(0, 0.15),
        initial = atanh(0.056),
        fixed = FALSE
      )
    )
  ) +

  f(COartime,
    model = "ar",
    order = 2,
    constr = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(1, 0.01),
        initial = log(1.118),
        fixed = FALSE
      ),
      pacf1 = list(
        prior = "normal",
        param = c(0, 0.15),
        initial = atanh(0.047),
        fixed = FALSE
      ),
      pacf2 = list(
        prior = "normal",
        param = c(0, 0.15),
        initial = atanh(0.056),
        fixed = FALSE
      )
    )
  ) +

  # =====================================================
  # Campo espaço-temporal SPDE
  # =====================================================

  f(s, 
    model=spde,group=s.group,
    control.group=list(model='ar1'),
    hyper = list(theta1 = list(initial = -4.148,
                               fixed = F),
                 theta2 = list(initial = -0.1267, 
                               fixed = F)))


resm2r <- inla(model2r, 
              family='poisson',
              data=inla.stack.data(stk), 
              E=exposure,
              control.predictor=list(A=inla.stack.A(stk), 
                                     compute=T),
              verbose=T, 
              control.inla = list(stupid.search=F,
                                  int.strategy = "eb"),num.threads=12)

save.image("modeloextremerain200interceptMulti.Rdata")							 
