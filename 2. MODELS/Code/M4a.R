library(spatstat)
library(sp)
library(gstat)
library(sf)
library(dplyr)
library(inlabru)
library(MatrixModels)

#setwd("Z:/home")

library(remotes)
library(INLA)
library(fmesher)
library(ggplot2)

inla.setOption(smtp = "taucs",  inla.mode = "classic", num.threads = "10:6")

num.threads <- "10:6"

set.seed(23)


load(file = "DATA/Data/FINAL ALL.RData")

all <- subset(all, is.na(`Site Name`) | PM25 >= 0)
sitedata <- subset(all, !is.na(`Site Name`))

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
model.covariates <- c("Time", "PM25PCM", "PM25MO_2", "background", "major_dens.scaled", 
                      "NDVI_New")

#Stacks
coordinates.allyear <- as.matrix(cbind(sitedata[,"Easting"], sitedata[,"Northing"]))

A_est <- inla.spde.make.A(mesh = mesh.s, loc = coordinates.allyear, group = sitedata$Time, n.group = 72)
s_index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde, n.group = 72)

A_SVC_est <- inla.spde.make.A(mesh = mesh.s, loc = coordinates.allyear,
                              weights = sitedata$PM25PCM)

svc_index <- inla.spde.make.index("spatial.varying.field", spde$n.spde)

names(sitedata)
stack_est <- inla.stack(data = list(logPM25 = sitedata$logPM25),
                        A = list(A_est, A_SVC_est, 1),
                        effects = list(s_index, svc_index, list(list = data.frame(Intercept = 1), 
                                                                sitedata[,model.covariates])), tag = "est")

i_month <- 72

A_pred <- inla.spde.make.A(mesh = mesh.s, loc = cbind(all$Easting, all$Northing), group = all$Time, n.group = 72)
A_SVC_pred <- inla.spde.make.A(mesh = mesh.s, loc = cbind(all$Easting, all$Northing),
                               weights = all$PM25PCM)

stack_pred <- inla.stack(data = list(logPM25 = NA), 
                         A = list(A_pred, A_SVC_pred, 1),
                         effects = list(s_index, svc_index, list(list = data.frame(Intercept = 1),
                                                                 all[,model.covariates])), tag="pred")

stack <- inla.stack(stack_est, stack_pred)



i.group <- 72


f.iid <- logPM25 ~ -1 + Intercept + background +  
  f(Time, PM25MO_2, model = "ar1", hyper = hyper.ar13) +
  f(spatial.field, model = spde, group = spatial.field.group, control.group = list(model="ar1", hyper = hyper.ar1.rho)) +
  f(spatial.varying.field, model = spde)



out <- inla(f.iid,
            family="gaussian",
            data=inla.stack.data(stack),
            control.predictor=list(A=inla.stack.A(stack), compute = T),
            verbose = F, num.threads = num.threads,
            control.inla = list(int.strategy = "eb"),
            control.compute = list(dic=TRUE, cpo=TRUE, waic=T, 
                                   mlik=T, return.marginals=T, config=T, 
                                   openmp.strategy="default", smtp="taucs"),            
            control.mode = list(theta = theta.ini, restart = TRUE)
            
)

save(out, file = "MODELS/MODEL/Outputs/M4a.RData")

load(file = "MODELS/MODEL/Outputs/M4a.RData")
summary(out)


index.est <- inla.stack.index(stack = stack, "est")$data
sitedata$pm25.estimated <- out$summary.fitted.values[index.est, "mean"]
cor(na.omit(sitedata[,c("logPM25","pm25.estimated")]))[1,2]

-mean(na.omit(out$cpo$cpo))

res <- sitedata$logPM25 - sitedata$pm25.estimated
sum(res^2) + sum(abs(res))


index.pred <- inla.stack.index(stack = stack, "pred")$data
all$pm25.predicted <- out$summary.fitted.values[index.pred, "mean"]

index.est <- inla.stack.index(stack = stack, "est")$data
sitedata$pm25.fitted <- out$summary.fitted.values[index.est, "mean"]

index.pred <- inla.stack.index(stack = stack, "pred")$data
all$pm25.predicted <- out$summary.fitted.values[index.pred, "mean"]
lp_summary <- out$summary.fitted.values[index.pred, "mean"]
lp_mean <- unlist(lapply(lp_summary, function(x) exp(x)))
all$pm25.predicted.exp <- lp_mean



###functions
RMSE <- function(z, zhat){
  z <- as.matrix(z)
  zhat <- as.matrix(zhat)
  res <- c(z - zhat)
  RMSE_val <- sqrt(mean(res^2,na.rm=T)) 
  return(RMSE_val)  
}

RMAE <- function(z, zhat){
  z <- as.matrix(z)
  zhat <- as.matrix(zhat)
  res <- c(z - zhat)
  RMAE_val <- sqrt(mean(abs(res) ,na.rm=T)) 
  return(RMAE_val)  
}

bias <- function(z, zhat){
  z <- as.matrix(z)
  zhat <- as.matrix(zhat)
  res <- c(zhat - z)
  bias_val <- mean(res) 
  return(bias_val)  
}

r2 <- function(z, zhat){
  z <- as.matrix(z)
  zhat <- as.matrix(zhat)
  res <- c(z - zhat)
  R2 <- cor(z,zhat)^2
  return(R2)
}

PMCC <- function(z, zhat){
  z <- as.matrix(z)
  zhat <- as.matrix(zhat)
  res <- c(z - zhat)
  gof <- sum(res^2)
  penalty <-  sum(abs(res))
  PMCC <- gof + penalty
  return(PMCC)
}

###Year CV
# 72 locs
set.seed(23)

coordinates.allyear <- as.matrix(cbind(sitedata[,"Easting"], sitedata[,"Northing"]))
sitedata$siteid <- as.numeric(as.factor(paste(sitedata$`Site Name`, sitedata$Network, sep = "_")))
sort(unique(sitedata$siteid))
n.sites <- length(unique(sitedata$siteid))

set.seed(23)

coords <- unique(sitedata[,c("Site Name", "Network", "siteid", "Easting", "Northing")])

r2_train_spde_a = r2_val_spde_a = rmse_train_spde_a = rmse_val_spde_a = c()
r2_train_spde_a_exp = r2_val_spde_a_exp = rmse_train_spde_a_exp = rmse_val_spde_a_exp = c()
r2_train_a = r2_val_a = se_train_a = se_val_a = ae_train_a = ae_val_a = var_train_a = var_val_a = c()
r2_train_a.exp = r2_val_a.exp = se_train_a.exp = se_val_a.exp = ae_train_a.exp = ae_val_a.exp = var_train_a.exp = var_val_a.exp =c()
cov <- c()

s_index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde, n.group =  72)

f.m1 <- logPM25 ~ -1 + Intercept + background +  
  f(Time, PM25MO_2, model = "ar1", hyper = hyper.ar13) +
  f(spatial.field, model = spde, group = spatial.field.group, control.group = list(model="ar1", hyper = hyper.ar1.rho)) +
  f(spatial.varying.field, model = spde)

outcome <- "logPM25"
results_list <- list()

out <- list()
for(j in 1:6){
  sitedata.train <- subset(sitedata, !(Year == j + 2013))
  sitedata.test <- subset(sitedata, Year == j + 2013)
  
  
  coords.train <- as.matrix(cbind(sitedata.train[,"Easting"], sitedata.train[,"Northing"]))
  coords.test <- as.matrix(cbind(sitedata.test[,"Easting"], sitedata.test[,"Northing"]))
  
  A_train <- inla.spde.make.A(mesh = mesh.s, loc = coords.train, group = sitedata.train$Time, n.group =  72)
  A_SVC_train <- inla.spde.make.A(mesh = mesh.s, loc = coords.train,
                                  weights = sitedata.train$PM25PCM)
  
  stack_train <- inla.stack(data = list(logPM25 = sitedata.train$logPM25),
                            A = list(A_train, A_SVC_train, 1),
                            effects = list(s_index, svc_index, list(list = data.frame(Intercept = 1), 
                                                                    sitedata.train[,model.covariates])),
                            tag = "train")
  
  
  i_month <-  72
  
  A_test <- inla.spde.make.A(mesh = mesh.s, loc = coords.test, group = sitedata.test$Time, n.group =  72)
  A_SVC_test <- inla.spde.make.A(mesh = mesh.s, loc = coords.test,
                                 weights = sitedata.test$PM25PCM)
  
  stack_test <- inla.stack(data = list(logPM25 = NA), 
                           A = list(A_test, A_SVC_test, 1),
                           effects = list(s_index, svc_index, list(list = data.frame(Intercept = 1),
                                                                   sitedata.test[,model.covariates])),
                           tag="test")
  
  
  
  stack <- inla.stack(stack_train, stack_test)
  i.group <-  72
  
  
  out[[j]] <- inla(f.m1,
                   family="gaussian",
                   data=inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack), compute=T),
                   verbose=F, num.threads = "8:4",
                   control.inla = list(int.strategy = "eb"),
                   control.compute = list(dic=TRUE, cpo=TRUE, waic=T, 
                                          mlik=T, return.marginals=T, config=T, 
                                          openmp.strategy="default", smtp="taucs"),
                   control.mode = list(theta = theta.ini, restart = TRUE)
                   
  )
  
  summary(out[[j]])
  index_train <- inla.stack.index(stack,"train")$data
  index_val <- inla.stack.index(stack,"test")$data
  
  results.train <- out[[j]]$summary.fitted$mean[index_train]
  results.val <- out[[j]]$summary.fitted$mean[index_val]
  
  
  cov <- out[[j]]$summary.fitted.values[index_val,] %>%
    bind_cols(obs = sitedata.test$logPM25) %>%
    mutate(obs_in_pi = obs >= `0.025quant` & obs <= `0.975quant`) %>%
    summarise(picp = sum(obs_in_pi) / n())
  
  # DIC and WAIC
  dic <- out[[j]]$waic$waic
  waic <- out[[j]]$dic$dic
  
  # Extract the fitted values
  M_fit <- out[[j]]$summary.fitted.values[,"mean"]
  
  # Squared Error
  se_train <- (sitedata.train$logPM25 - results.train)^2
  se_val <- (sitedata.test$logPM25 - results.val)^2
  
  # R2: Pseudo_r2 function defined above
  r2_train <- cor(sitedata.train$logPM25, results.train)^2
  r2_val <- cor(sitedata.test$logPM25, results.val)^2
  
  #
  rmse_train <- RMSE(sitedata.train$logPM25, results.train)
  rmse_val <- RMSE(sitedata.test$logPM25, results.val)
  rmae_train <- RMAE(sitedata.train$logPM25, results.train)
  rmae_val <- RMAE(sitedata.test$logPM25, results.val)
  bias_train <- bias(sitedata.train$logPM25, results.train)
  bias_val <- bias(sitedata.test$logPM25, results.val)
  r2_train <- r2(sitedata.train$logPM25, results.train)
  r2_val <- r2(sitedata.test$logPM25, results.val)
  pmcc_train <- PMCC(sitedata.train$logPM25, results.train)
  pmcc_val <- PMCC(sitedata.test$logPM25, results.val)
  
  
  # Store the results
  results_list[[j]] = list("dic" = dic, # 1 - 5
                           "waic" = waic, # 6 - 10
                           "se_val" = se_val, # 11
                           "se_train" = se_train, # 16
                           "r2_val" = r2_val, # 21
                           "r2_train" = r2_train, # 26
                           "rmse_train" = rmse_train, # 31
                           "rmse_val" = rmse_val, # 36
                           "rmae_train" = rmae_train, # 41
                           "rmae_val" = rmae_val, # 46
                           "bias_train" = bias_train, # 51
                           "bias_val" = bias_val, # 56
                           "pmcc_train" = pmcc_train, # 61
                           "pmcc_val" =  pmcc_val, # 66
                           "cov" = cov) # 71
  
  
}    

results <- do.call(rbind, results_list)
save(results, file = "MODELS/MODEL/Outputs/M4a_TimeCV.RData")

load(file = "MODELS/MODEL/Outputs/M4a_TimeCV.RData")



results

# R2
mean(unlist(results[c(25:30)])) #5

# RMSE
mean(unlist(results[c(43:48)])) #8

# bias
mean(unlist(results[c(67:72)])) #12

# PMCC
mean(unlist(results[c(79:84)])) #14

# Cov
mean(unlist(results[c(85:90)])) #15






## Spatial Block CV
sitedata$siteid <- as.numeric(as.factor(paste(sitedata$`Site Name`, sitedata$Network)))
sort(unique(sitedata$siteid))

coords <- unique(sitedata[,c("Site Name", "siteid", "Easting", "Northing")])

a <- quantile(all$Northing, c(0:2/2))[2]
b1 <- quantile(all$Easting, c(0:3/3))[2]
b2 <- quantile(all$Easting, c(0:3/3))[3]


plot(coords[,c(3,4)])
abline(h = a)
abline(v = b1)
abline(v = b2)

sitedata.test.sets <- list()
sitedata.train.sets <- list()

sitedata.test.sets[[1]] <- subset(sitedata, Easting < b1 & Northing < a)
sitedata.train.sets[[1]] <- subset(sitedata, !(Easting < b1 & Northing < a))

sitedata.test.sets[[2]] <- subset(sitedata, Easting < b1 & Northing > a)
sitedata.train.sets[[2]] <- subset(sitedata, !(Easting < b1 & Northing > a))

sitedata.test.sets[[3]] <- subset(sitedata, Easting > b1 & Easting < b2 & Northing < a)
sitedata.train.sets[[3]] <- subset(sitedata, !(Easting > b1 & Easting < b2 & Northing < a))

sitedata.test.sets[[4]] <- subset(sitedata, Easting > b1 & Easting < b2 & Northing > a)
sitedata.train.sets[[4]] <- subset(sitedata, !(Easting > b1 & Easting < b2 & Northing > a))

sitedata.test.sets[[5]] <- subset(sitedata, Easting > b2 & Northing < a)
sitedata.train.sets[[5]] <- subset(sitedata, !(Easting > b2 & Northing < a))

sitedata.test.sets[[6]] <- subset(sitedata, Easting > b2 & Northing > a)
sitedata.train.sets[[6]] <- subset(sitedata, !(Easting > b2 & Northing > a))


r2_train_spde_a = r2_val_spde_a = rmse_train_spde_a = rmse_val_spde_a = c()
r2_train_spde_a_exp = r2_val_spde_a_exp = rmse_train_spde_a_exp = rmse_val_spde_a_exp = c()
r2_train_a = r2_val_a = se_train_a = se_val_a = ae_train_a = ae_val_a = var_train_a = var_val_a = c()
r2_train_a.exp = r2_val_a.exp = se_train_a.exp = se_val_a.exp = ae_train_a.exp = ae_val_a.exp = var_train_a.exp = var_val_a.exp =c()
cov <- c()

s_index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde, n.group =  72)


f.m1 <- logPM25 ~ -1 + Intercept + background +  
  f(Time, PM25MO_2, model = "ar1", hyper = hyper.ar13) +
  f(spatial.field, model = spde, group = spatial.field.group, control.group = list(model="ar1", hyper = hyper.ar1.rho)) +
  f(spatial.varying.field, model = spde)

outcome <- "logPM25"
results_list <- list()

out <- list()
for(j in 1:6){
  sitedata.train <- subset(sitedata, siteid %in% sitedata.train.sets[[j]]$siteid)
  sitedata.test <- subset(sitedata, siteid %in% sitedata.test.sets[[j]]$siteid)
  
  coords.train <- as.matrix(cbind(sitedata.train[,"Easting"], sitedata.train[,"Northing"]))
  coords.test <- as.matrix(cbind(sitedata.test[,"Easting"], sitedata.test[,"Northing"]))
  
  A_train <- inla.spde.make.A(mesh = mesh.s, loc = coords.train, group = sitedata.train$Time, n.group =  72)
  A_SVC_train <- inla.spde.make.A(mesh = mesh.s, loc = coords.train,
                                  
                                  weights = sitedata.train$PM25PCM)
  
  stack_train <- inla.stack(data = list(logPM25 = sitedata.train$logPM25),
                            A = list(A_train, A_SVC_train, 1),
                            effects = list(s_index, svc_index, list(list = data.frame(Intercept = 1), 
                                                                    sitedata.train[,model.covariates])),
                            tag = "train")
  
  
  i_month <-  72
  
  A_test <- inla.spde.make.A(mesh = mesh.s, loc = coords.test, group = sitedata.test$Time, n.group =  72)
  A_SVC_test <- inla.spde.make.A(mesh = mesh.s, loc = coords.test,
                                 
                                 weights = sitedata.test$PM25PCM)
  
  stack_test <- inla.stack(data = list(logPM25 = NA), 
                           A = list(A_test, A_SVC_test, 1),
                           effects = list(s_index, svc_index, list(list = data.frame(Intercept = 1),
                                                                   sitedata.test[,model.covariates])),
                           tag="test")
  
  
  
  
  stack <- inla.stack(stack_train, stack_test)
  i.group <-  72
  
  
  out[[j]] <- inla(f.m1,
                   family="gaussian",
                   data=inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack), compute=T),
                   verbose=F, num.threads = "8:4",
                   control.inla = list(int.strategy = "eb"),
                   control.compute = list(dic=TRUE, cpo=TRUE, waic=T, 
                                          mlik=T, return.marginals=T, config=T, 
                                          openmp.strategy="default", smtp="taucs"),
                   control.mode = list(theta = theta.ini, restart = TRUE)
                   
  )
  summary(out[[j]])
  index_train <- inla.stack.index(stack,"train")$data
  index_val <- inla.stack.index(stack,"test")$data
  
  results.train <- out[[j]]$summary.fitted$mean[index_train]
  results.val <- out[[j]]$summary.fitted$mean[index_val]
  
  
  cov <- out[[j]]$summary.fitted.values[index_val,] %>%
    bind_cols(obs = sitedata.test$logPM25) %>%
    mutate(obs_in_pi = obs >= `0.025quant` & obs <= `0.975quant`) %>%
    summarise(picp = sum(obs_in_pi) / n())
  
  # DIC and WAIC
  dic <- out[[j]]$waic$waic
  waic <- out[[j]]$dic$dic
  
  # Extract the fitted values
  M_fit <- out[[j]]$summary.fitted.values[,"mean"]
  
  # Squared Error
  se_train <- (sitedata.train$logPM25 - results.train)^2
  se_val <- (sitedata.test$logPM25 - results.val)^2
  
  # R2: Pseudo_r2 function defined above
  r2_train <- cor(sitedata.train$logPM25, results.train)^2
  r2_val <- cor(sitedata.test$logPM25, results.val)^2
  
  #
  rmse_train <- RMSE(sitedata.train$logPM25, results.train)
  rmse_val <- RMSE(sitedata.test$logPM25, results.val)
  rmae_train <- RMAE(sitedata.train$logPM25, results.train)
  rmae_val <- RMAE(sitedata.test$logPM25, results.val)
  bias_train <- bias(sitedata.train$logPM25, results.train)
  bias_val <- bias(sitedata.test$logPM25, results.val)
  r2_train <- r2(sitedata.train$logPM25, results.train)
  r2_val <- r2(sitedata.test$logPM25, results.val)
  pmcc_train <- PMCC(sitedata.train$logPM25, results.train)
  pmcc_val <- PMCC(sitedata.test$logPM25, results.val)
  
  
  # Store the results
  results_list[[j]] = list("dic" = dic,
                           "waic" = waic, 
                           "se_val" = se_val,
                           "se_train" = se_train,
                           "r2_val" = r2_val,
                           "r2_train" = r2_train,
                           "rmse_train" = rmse_train, 
                           "rmse_val" = rmse_val,
                           "rmae_train" = rmae_train,
                           "rmae_val" = rmae_val,
                           "bias_train" = bias_train,
                           "bias_val" = bias_val,
                           "pmcc_train" = pmcc_train,
                           "pmcc_val" =  pmcc_val,
                           "cov" = cov)
  
  
}    

results <- do.call(rbind, results_list)
save(results, file = "MODELS/MODEL/Outputs/M4a_SpatialBlockCV.RData")

#load(file = "MODELS/Outputs/M3_New_NDVI_SpatialBlockCV.RData")

results

# R2
mean(unlist(results[c(25:30)])) #5

# RMSE
mean(unlist(results[c(43:48)])) #8

# bias
mean(unlist(results[c(67:72)])) #12

# PMCC
mean(unlist(results[c(79:84)])) #14

# Cov
mean(unlist(results[c(85:90)])) #15











