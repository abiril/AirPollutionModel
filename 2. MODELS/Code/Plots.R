library(spatstat)
library(sp)
library(gstat)
library(sf)
library(dplyr)
library(inlabru)
library(MatrixModels)

setwd("Z:/home")

library(remotes)
library(INLA)
library(fmesher)
library(ggplot2)

inla.setOption(smtp = "taucs",  inla.mode = "classic", num.threads = "8:8")

num.threads <- "8:8"

set.seed(23)


load(file = "DATA/Data/FINAL ALL.RData")

all <- subset(all, is.na(`Site Name`) | PM25 >= 0)
sitedata <- subset(all,  PM25 >= 0)

# Make mesh

bound <- inla.nonconvex.hull(as.matrix(all[,c("Easting", "Northing")]),
                             convex = 5)

mesh.s <- fmesher::fm_mesh_2d_inla(
  locs = cbind(sitedata$Easting, sitedata$Northing),
  boundary = bound,
  offset = c(0.1,5), max.edge = c(2.5,4))


# Parameters and priors

t.max <- 72
mesh.t <- inla.mesh.1d(1:t.max)

spde <- inla.spde2.pcmatern(
  mesh = mesh.s,
  prior.range = c(20, 0.01),
  prior.sigma = c(0.1, 0.1),
  alpha = 3/2
)

hyper.iid.fix <- list(prec = list(initial = -2*log(0.1)))
hyper2 <- list(prec = list(prior="pc.prec", param=c(0.5,0.5)))
hyper.ar1.rho <-list(rho = list(prior = "pc.cor1", param=c(0.95,0.5)))
hyper22 <- list(prec = list(initial = -2*log(5)))
hyper.ar1.rho2 <- list(rho = list(prior = "pc.cor1", param=c(0.5,0.5))) # P(rho > 0.95) = 0.5

hyper.ar13 <- list(prec = list(param = c(5,5)))


# Initial Values

theta.ini <- c(
  log(50), # likelihood log precision
  log(10), # tvc precision
  log((1 + 0.85)/(1 - 0.85)), # log trans rho tvc
  log(20), # log range
  log(0.35), # log stdev
  log((1 + 0.95)/(1 - 0.95)), # log trans rho 1
  log(10), # log range for SVC - trying to push this down
  log(0.03)) # log stdev for SVC ) # log trans rho for SVC


#INLA Set-Up
model.covariates <- c("Time", "PM25PCM", "PM25MO_2", "background", "all_dens.scaled", 
                      "NDVI_New")

#Stacks
coordinates.allyear <- as.matrix(cbind(sitedata[,"Easting"], sitedata[,"Northing"]))

A_est <- inla.spde.make.A(mesh = mesh.s, loc = coordinates.allyear, group = sitedata$Time, n.group = 72)
s_index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde, n.group = 72)

A_SVC_est <- inla.spde.make.A(mesh = mesh.s, loc = coordinates.allyear,
                              weights = sitedata$PM25MO_2)

svc_index <- inla.spde.make.index("spatial.varying.field", spde$n.spde)

names(sitedata)
stack_est <- inla.stack(data = list(logPM25 = sitedata$logPM25),
                        A = list(A_est, A_SVC_est, 1),
                        effects = list(s_index, svc_index, list(list = data.frame(Intercept = 1), 
                                                                sitedata[,model.covariates])), tag = "est")

i_month <- 72

A_pred <- inla.spde.make.A(mesh = mesh.s, loc = cbind(all$Easting, all$Northing), group = all$Time, n.group = 72)
A_SVC_pred <- inla.spde.make.A(mesh = mesh.s, loc = cbind(all$Easting, all$Northing),
                               weights = all$PM25MO_2)

stack_pred <- inla.stack(data = list(logPM25 = NA), 
                         A = list(A_pred, A_SVC_pred, 1),
                         effects = list(s_index, svc_index, list(list = data.frame(Intercept = 1),
                                                                 all[,model.covariates])), tag="pred")

stack <- inla.stack(stack_est, stack_pred)

i.group <- 72


f.iid <- logPM25 ~ -1 + Intercept + background + PM25PCM +  
  f(Time, NDVI_New, model = "ar1", hyper = hyper.ar13) +
  f(spatial.field, model = spde, group = spatial.field.group, control.group = list(model="ar1", hyper = hyper.ar1.rho)) +
  f(spatial.varying.field, model = spde)

load(file = "MODELS/MODEL/Outputs/M4b_Laplace.RData")




library(ggopenair)

setwd("Z:/home")

index.y <- inla.stack.index(stack = stack, "est")$data
index.pred <- inla.stack.index(stack = stack, "pred")$data

all$pm25.predicted <- out$summary.fitted.values[index.pred, "mean"]
all$pm25.predicted.sd <- out$summary.fitted.values[index.pred, "sd"]

lp_summary <- out$summary.fitted.values[index.pred, "mean"]
lp_mean <- unlist(lapply(lp_summary, function(x) exp(x)))
all$pm25.predicted.exp <- lp_mean

lp_summary <- out$summary.fitted.values[index.y, "mean"]
lp_mean <- unlist(lapply(lp_summary, function(x) exp(x)))
sitedata$pm25.predicted.exp <- lp_mean

lp_summary <- out$summary.fitted.values[index.pred, "sd"]
lp_mean <- unlist(lapply(lp_summary, function(x) exp(x)))
all$pm25.predicted.sd.exp <- lp_mean

all$x <- all$Easting
all$y <- all$Northing

all.grid <- subset(all, is.na(siteid))

all.sf <- st_as_sf(all.grid, coords = c("Easting","Northing"), 
                   crs = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.999601272 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")



#month.i <- c(13, 20, 49)
month.i <- c(13, 22, 42)
month.i <- c(13, 34, 65)
month.i <- c(12, 33, 65)
prob <- 10

# Pred PM25
all.month <- subset(all.sf, Time %in% month.i)

scale.col <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf",
               "#fee090","#fdae61","#f46d43","#d73027","#a50026")

quant.all <- as.numeric(unlist(quantile(all.month$pm25.predicted.exp, seq(0, 1, 0.09090909))))
quant.all

quant.all <- c(4.5, 7.3, 7.9, 8.4, 8.9, 9.4, 9.9, 10.6, 11.5, 13.0, 16.3, 40.8) #c(13, 22)

quant.all <- c(4.5, 7, 8, 8.5, 9, 9.5, 10, 10.5, 11.5, 13.0, 16.5, 41) #c(13, 24)
quant.all <- c(4.9, 6, 7, 8, 8.5, 9, 9.5, 10, 10.5, 12, 15, 30) #c(13, 22, 56)
quant.all <- c(6.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 15, 30) #c(13, 22, 65)


quant.all <- as.numeric(unlist(quantile(all.month$pm25.predicted.exp, seq(0, 1, 0.08333333))))
quant.all
quant.all <- c(5, 8, 8.5, 9, 9.25, 9.5, 9.75, 10, 10.33, 10.66, 11, 11.5, 12, 14, 23) #c(12, 33, 65)

#quant.all <- c(5, 8, 8.5, 9, 9.25, 9.5, 9.75, 10, 10.25, 10.5, 10.75, 11, 11.5, 12, 14, 23) #c(12, 33, 65)

#quant.all <- c(5, 9, 9.5, 9.75, 10, 10.25, 10.5, 10.75, 11, 11.5, 12, 20)


p.predPM25 <- ggplot() +
  geom_raster(data = all.month, aes(x = x, y = y, fill = pm25.predicted.exp)) +
  scale_fill_stepsn(
    colours = scale.col,
    breaks = quant.all,
    limits = c(quant.all[1], quant.all[15]),
    values = scales::rescale(quant.all)
  ) +
  facet_wrap(~ Time, ncol = 1, nrow = length(month.i)) +
  theme(panel.spacing = unit(0, "points"),
        strip.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_blank() , 
        strip.text.x = element_blank(),
        legend.position = "right",
        plot.title = element_text(size=30),
        legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 14)) +
  coord_equal() +
  labs(fill = "PM2.5") +
  ggtitle("Predicted PM2.5")


# AQR

load(file = "DATA/data/20-Met Office Reanalysis PM25 Cubic.RData")

AQR.month <- subset(PM25.DF, Time %in% (month.i + 48))
range(AQR.month$PM25)

p.aqrPM25 <- ggplot() +
  geom_raster(data = AQR.month, aes(x = x, y = y, fill = PM25)) +
  scale_fill_stepsn(
    colours = scale.col,
    breaks = quant.all,
    limits = c(quant.all[1], quant.all[15]),
    values = scales::rescale(quant.all)
  ) +  
  facet_wrap(~ Time, ncol = 1, nrow = length(month.i)) +
  theme(panel.spacing = unit(0, "points"),
        strip.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_blank() , 
        strip.text.x = element_blank(),
        legend.position = "right",
        plot.title = element_text(size=30),
        legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 14)) +
  coord_equal() +
  labs(fill = "PM2.5")+
  ggtitle("AQR Model PM2.5")



# SD Pred 

breaks <- as.numeric(c(quantile(na.omit(all.sf$pm25.predicted.sd.exp), probs = seq(0, 1, by = .1))))

p.sdPM25 <- ggplot() +
  geom_raster(data = all.month, aes(x = x, y = y, fill = pm25.predicted.sd.exp)) +
  #scale_fill_viridis(option="magma", na.value = "transparent", breaks = breaks, limits=range(breaks)) +
  scale_fill_viridis(option="C", na.value = "transparent") +
  facet_wrap(~ Time, ncol = 1, nrow = length(month.i)) +
  theme(panel.spacing = unit(0, "points"),
        strip.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_blank() , 
        strip.text.x = element_blank(),
        legend.position = "right",
        plot.title = element_text(size=30),
        legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 14)) +
  coord_equal() +
  labs(fill = "SD") +
  ggtitle("Predicted PM2.5 SD")


# Pred Prob PM25 > 10
scale.col <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf",
               "#fee090","#fdae61","#f46d43","#d73027","#a50026")

all.month$predProb = pnorm((all.month$pm25.predicted.exp - prob)/sd(all.month$pm25.predicted.exp))
all.month$excProb = cut2(all.month$predProb, g=3, cuts = c(0.2,0.8))
quant.all <- as.numeric(unlist(quantile(all.month$predProb, seq(0, 1, 0.09090909))))

p.predProbPM25 <- ggplot() +
  geom_raster(data = all.month, aes(x = x, y = y, fill = excProb)) +
  scale_fill_viridis_d(option = "viridis", direction = 1, begin = 0.2, end = 1,
                       labels = c("[0, 0.2)", "[0.2, 0.8)", "[0.8,1]")) +
  facet_wrap(~ Time, ncol = 1, nrow = length(month.i)) +
  theme(panel.spacing = unit(0, "points"),
        strip.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_blank() , 
        strip.text.x = element_blank(),
        legend.position = "right",
        plot.title = element_text(size=30),
        legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 14)) +
  coord_equal() +
  labs(fill = "Exc Prob") +
  ggtitle("Exceedance Probability PM2.5 > 10")



# Spatial Field
projgrid <- inla.mesh.projector(mesh.s, xlim = range(all.grid$Easting), 
                                ylim = range(all.grid$Northing), dims = c(62,49))

library(MatrixModels)
month.spatial.df <- list()
for(i in 1:length(month.i)){
  xmean <- inla.mesh.project(projgrid, out$summary.random$spatial.field$mean[s_index$spatial.field.group == month.i[i]])
  
  month.spatial.field <- list(x = projgrid$x, y = projgrid$y, z = xmean)
  
  month.spatial.field.r <- terra::rast(month.spatial.field, warn = TRUE)
  
  month.spatial.df[[i]] <- as.data.frame(month.spatial.field.r, xy = TRUE)
  names(month.spatial.df[[i]])[3] <- "spatial.field"
  month.spatial.df[[i]]$Time <- month.i[i]
  
}

spatial.df <- do.call(rbind, month.spatial.df)

library(paletteer)

scale.col <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf",
               "#fee090","#fdae61","#f46d43","#d73027","#a50026")


p.spatialField <- ggplot() +
  geom_raster(spatial.df, mapping = aes(x = x, y = y, fill = spatial.field)) +
  scale_fill_paletteer_c("grDevices::Purple-Yellow") +
  facet_wrap(~ Time, ncol = 1, nrow = length(month.i)) +
  theme(panel.spacing = unit(0, "points"),
        strip.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_blank(), 
        strip.text.x = element_blank(),
        legend.position = "right",
        plot.title = element_text(size=30),
        legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 14)) +
  coord_equal() +
  labs(fill = "Sp. Field") +
  ggtitle("Spatial Field for SPDE")

########

ggarrange(p.aqrPM25, p.predPM25, p.sdPM25, p.spatialField,  p.predProbPM25,
          nrow = 1, ncol = 5,
          align = "hv")

ggsave("MODELS/MODEL/Plots/M4b 5 Paper Plots 3 Months.png",
       width = 75, height = 60, units = "cm")

