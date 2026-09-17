library(INLA)
library(sf)
library(fmesher)
library(sp)
setwd("/dados/extremerain")

load("modeloextremerain.Rdata")

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



stk <- inla.stack(data=list(y=agg.dat$Freq, exposure=e0),
                  A=list(A.st, 1),
                  effects=list(idx,
                               list(Intercept=rep(1,
                                           nrow(agg.dat)),
                                    rwtime=seastime,
                                    seastime=seastime,
                                    artime=seastime)))




model2 <- y ~ -1 + Intercept+
  f(rwtime,
    model = "rw2",
    constr = TRUE,
    scale.model = TRUE,
    hyper = list(
      prec = list(
        prior = "pc.prec",
        param = c(10, 0.1),
        initial = log(118.253)
      )
    )
  ) +

  f(seastime,
    model = "seasonal",
    season.length = 12,
    constr = TRUE,
    hyper = list(
      prec = list(
        initial = log(278.591),
        fixed = FALSE
      )
    )
  ) +

  f(artime,
    model = "ar",
    order = 2,
    constr = TRUE,
    hyper = list(
      prec  = list(
        initial = log(1.118),
        fixed = FALSE
      ),
      pacf1 = list(
        initial = atanh(0.047),
        fixed = FALSE
      ),
      pacf2 = list(
        initial = atanh(0.056),
        fixed = FALSE
      )
    )
  )+
  f(s, 
    model=spde,group=s.group,
    control.group=list(model='ar1'),
    hyper = list(theta1 = list(initial = -4.148,
                               fixed = F),
                 theta2 = list(initial = -0.1267, 
                               fixed = F)))

resm2 <- inla(model2, 
              family='poisson',
              data=inla.stack.data(stk), 
              E=exposure,
              control.predictor=list(A=inla.stack.A(stk), 
                                     compute=T),
              verbose=T, 
              control.inla = list(stupid.search=F,
                                  int.strategy = "eb"))

save.image("modeloextremerain100intercept.Rdata")							 
