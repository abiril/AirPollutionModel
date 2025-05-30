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


#INLA Set-Up
model.covariates <- c("Time", "PM25PCM", "PM25MO_2", "background", "major_dens.scaled", "NDVI_New")

cor(sitedata[,c("logPM25", model.covariates)], use = "pairwise.complete.obs")
cor(all[,model.covariates], use = "pairwise.complete.obs")

coordinates.allyear <- as.matrix(cbind(sitedata[,"Easting"], sitedata[,"Northing"]))

est <- inla.spde.make.A(mesh = mesh.s, loc = coordinates.allyear, group = sitedata$Time, n.group = 72)
s_index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde, n.group = 72)

names(sitedata)
stack_est <- inla.stack(data = list(logPM25 = sitedata$logPM25),
                        A = list(est, 1),
                        effects = list(s_index, list(list = data.frame(Intercept = 1), 
                                                     sitedata[,model.covariates])), tag = "est")

i_month <- 72

pred <- inla.spde.make.A(mesh = mesh.s, loc = cbind(all$Easting, all$Northing), group = all$Time, n.group = 72)
stack_pred <- inla.stack(data = list(logPM25 = NA), 
                         A = list(pred, 1),
                         effects = list(s_index, list(list = data.frame(Intercept = 1),
                                                      all[,model.covariates])), tag="pred")

stack <- inla.stack(stack_est, stack_pred)

i.group <- 72


#################################
#################################
########## MODEL ################
#################################
#################################

# "PM25PCM", "PM25AQUM", "Temp", "background", "NDVI", "pop"

f.iid <- logPM25 ~ -1 + Intercept +  PM25PCM + PM25MO_2 + background +   
  f(spatial.field, model = spde, group = spatial.field.group, control.group = list(model="ar1", hyper = hyper.ar1.rho))

out <- inla(f.iid,
            family="gaussian",
            data=inla.stack.data(stack),
            control.predictor=list(A=inla.stack.A(stack), compute = T),
            verbose = F, num.threads = num.threads,
            control.inla = list(int.strategy = "eb"),
            control.compute = list(dic=TRUE, cpo=TRUE, waic=T, 
                                   mlik=T, return.marginals=T, config=T, 
                                   openmp.strategy="default", smtp="taucs")
)

save(out, file = "MODELS/MODEL/Outputs/M1a.RData")

load(file = "MODELS/MODEL/Outputs/M1a.RData")
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
all$pm25.predicted_exp <- lp_mean




#################################
#################################
####### Cross-Val Time ##########
########## M1 Both ##############
#################################
#################################

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

r2_train_spde = r2_val_spde = rmse_train_spde = rmse_val_spde = c()
r2_train_spde_exp = r2_val_spde_exp = rmse_train_spde_exp = rmse_val_spde_exp = c()

r2_train = r2_val = se_train = se_val = ae_train = ae_val = var_train = var_val = c()
r2_train_exp = r2_val_exp = se_train_exp = se_val_exp = ae_train_exp = ae_val_exp = var_train_exp = var_val_exp =c()
cov <- c()

s_index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde, n.group =  72)

f.m1 <- logPM25 ~ -1 + Intercept +  PM25MO_2 + PM25PCM + background +   
  f(spatial.field, model = spde, group = spatial.field.group, control.group = list(model="ar1", hyper = hyper.ar1.rho))

outcome <- "logPM25"
results_list <- list()

out <- list()
for(j in 1:6){
  sitedata.train <- subset(sitedata, !(Year == j + 2013))
  sitedata.test <- subset(sitedata, Year == j + 2013)
  
  
  coords.train <- as.matrix(cbind(sitedata.train[,"Easting"], sitedata.train[,"Northing"]))
  coords.test <- as.matrix(cbind(sitedata.test[,"Easting"], sitedata.test[,"Northing"]))
  
  train <- inla.spde.make.A(mesh = mesh.s, loc = coords.train, group = sitedata.train$Time, n.group =  72)
  
  stack_train <- inla.stack(data = list(logPM25 = sitedata.train$logPM25),
                            A = list(train, 1),
                            effects = list(s_index, list(list = data.frame(Intercept = 1), 
                                                         sitedata.train[,model.covariates])),
                            tag = "train")
  
  
  i_month <-  72
  
  test <- inla.spde.make.A(mesh = mesh.s, loc = coords.test, group = sitedata.test$Time, n.group =  72)
  
  stack_test <- inla.stack(data = list(logPM25 = NA), 
                           A = list(test, 1),
                           effects = list(s_index, list(list = data.frame(Intercept = 1),
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
                                          openmp.strategy="default", smtp="taucs")
  )
  summary(out[[j]])
  index_train <- inla.stack.index(stack,"train")$data
  index_val <- inla.stack.index(stack,"test")$data
  
  results.train <- out[[j]]$summary.fitted$mean[index_train]
  results.val <- out[[j]]$summary.fitted$mean[index_val]
  
  results.train_exp <- unlist(lapply(out[[j]]$summary.fitted$mean[index_train], function(x) exp(x)))
  results.val_exp <- unlist(lapply(out[[j]]$summary.fitted$mean[index_val], function(x) exp(x)))
  
  
  cov <- out[[j]]$summary.fitted.values[index_val,] %>%
    bind_cols(obs = sitedata.test$logPM25) %>%
    mutate(obs_in_pi = obs >= `0.025quant` & obs <= `0.975quant`) %>%
    summarise(picp = sum(obs_in_pi) / n())
  
  # DIC and WAIC
  dic <- out[[j]]$waic$waic
  waic <- out[[j]]$dic$dic
  
  # Extract the fitted values
  M_fit <- out[[j]]$summary.fitted.values[,"mean"]
  M_fit_exp <- unlist(lapply(out[[j]]$summary.fitted.values[,"mean"], function(x) exp(x)))
  
  # Squared Error
  se_train <- (sitedata.train$logPM25 - results.train)^2
  se_val <- (sitedata.test$logPM25 - results.val)^2
  
  se_train_exp <- (sitedata.train$PM25 - results.train)^2
  se_val_exp <- (sitedata.test$PM25 - results.val)^2
  
  # R2: Pseudo_r2 function defined above
  r2_train <- cor(sitedata.train$logPM25, results.train)^2
  r2_val <- cor(sitedata.test$logPM25, results.val)^2
  
  r2_train_exp <- cor(sitedata.train$PM25, results.train_exp)^2
  r2_val_exp <- cor(sitedata.test$PM25, results.val)^2
  
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
  
  rmse_train_exp <- RMSE(sitedata.train$PM25, results.train_exp)
  rmse_val_exp <- RMSE(sitedata.test$PM25, results.val_exp)
  rmae_train_exp <- RMAE(sitedata.train$PM25, results.train_exp)
  rmae_val_exp <- RMAE(sitedata.test$PM25, results.val_exp)
  bias_train_exp <- bias(sitedata.train$PM25, results.train_exp)
  bias_val_exp <- bias(sitedata.test$PM25, results.val_exp)
  r2_train_exp <- r2(sitedata.train$PM25, results.train_exp)
  r2_val_exp <- r2(sitedata.test$PM25, results.val_exp)
  pmcc_train_exp <- PMCC(sitedata.train$PM25, results.train_exp)
  pmcc_val_exp <- PMCC(sitedata.test$PM25, results.val_exp)
  
  
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
                           "cov" = cov, # 71
                           "se_val_exp" = se_val_exp, # 11
                           "se_train_exp" = se_train_exp, # 16
                           "r2_val_exp" = r2_val_exp, # 21
                           "r2_train_exp" = r2_train_exp, # 26
                           "rmse_train_exp" = rmse_train_exp, # 31
                           "rmse_val_exp" = rmse_val_exp, # 36
                           "rmae_train_exp" = rmae_train_exp, # 41
                           "rmae_val_exp" = rmae_val_exp, # 46
                           "bias_train_exp" = bias_train_exp, # 51
                           "bias_val_exp" = bias_val_exp, # 56
                           "pmcc_train_exp" = pmcc_train_exp, # 61
                           "pmcc_val_exp" =  pmcc_val_exp) # 66
  
  
  
}    

results <- do.call(rbind, results_list)
save(results, file = "MODELS/MODEL/Outputs/M1a_exp_TimeCV.RData")

load(file = "MODELS/MODEL/Outputs/M1a_exp_TimeCV.RData")

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

quantile(all$Easting, c(0:3/3))[2]
quantile(all$Northing, c(0:2/2))


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


r2_train_spde = r2_val_spde = rmse_train_spde = rmse_val_spde = c()
r2_train_spde_exp = r2_val_spde_exp = rmse_train_spde_exp = rmse_val_spde_exp = c()
r2_train = r2_val = se_train = se_val = ae_train = ae_val = var_train = var_val = c()
r2_train_exp = r2_val_exp = se_train_exp = se_val_exp = ae_train_exp = ae_val_exp = var_train_exp = var_val_exp =c()
cov <- c()

s_index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde, n.group =  72)


f.m1 <- logPM25 ~ -1 + Intercept +  PM25MO_2 + PM25PCM + background +   
  f(spatial.field, model = spde, group = spatial.field.group, control.group = list(model="ar1", hyper = hyper.ar1.rho))

outcome <- "logPM25"
results_list <- list()

out <- list()
for(j in 1:6){
  sitedata.train <- subset(sitedata, siteid %in% sitedata.train.sets[[j]]$siteid)
  sitedata.test <- subset(sitedata, siteid %in% sitedata.test.sets[[j]]$siteid)
  
  coords.train <- as.matrix(cbind(sitedata.train[,"Easting"], sitedata.train[,"Northing"]))
  coords.test <- as.matrix(cbind(sitedata.test[,"Easting"], sitedata.test[,"Northing"]))
  
  train <- inla.spde.make.A(mesh = mesh.s, loc = coords.train, group = sitedata.train$Time, n.group =  72)
  
  stack_train <- inla.stack(data = list(logPM25 = sitedata.train$logPM25),
                            A = list(train, 1),
                            effects = list(s_index, list(list = data.frame(Intercept = 1), 
                                                         sitedata.train[,model.covariates])),
                            tag = "train")
  
  
  i_month <-  72
  
  test <- inla.spde.make.A(mesh = mesh.s, loc = coords.test, group = sitedata.test$Time, n.group =  72)
  
  stack_test <- inla.stack(data = list(logPM25 = NA), 
                           A = list(test, 1),
                           effects = list(s_index, list(list = data.frame(Intercept = 1),
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
                                          openmp.strategy="default", smtp="taucs")
  )
  summary(out[[j]])
  index_train <- inla.stack.index(stack,"train")$data
  index_val <- inla.stack.index(stack,"test")$data
  
  results.train <- out[[j]]$summary.fitted$mean[index_train]
  results.val <- out[[j]]$summary.fitted$mean[index_val]
  
  results.train_exp <- unlist(lapply(out[[j]]$summary.fitted$mean[index_train], function(x) exp(x)))
  results.val_exp <- unlist(lapply(out[[j]]$summary.fitted$mean[index_val], function(x) exp(x)))
  
  
  cov <- out[[j]]$summary.fitted.values[index_val,] %>%
    bind_cols(obs = sitedata.test$logPM25) %>%
    mutate(obs_in_pi = obs >= `0.025quant` & obs <= `0.975quant`) %>%
    summarise(picp = sum(obs_in_pi) / n())
  
  
  # DIC and WAIC
  dic <- out[[j]]$waic$waic
  waic <- out[[j]]$dic$dic
  
  # Extract the fitted values
  M_fit <- out[[j]]$summary.fitted.values[,"mean"]
  M_fit_exp <- unlist(lapply(out[[j]]$summary.fitted.values[,"mean"], function(x) exp(x)))
  
  # Squared Error
  se_train <- (sitedata.train$logPM25 - results.train)^2
  se_val <- (sitedata.test$logPM25 - results.val)^2
  
  se_train_exp <- (sitedata.train$PM25 - results.train)^2
  se_val_exp <- (sitedata.test$PM25 - results.val)^2
  
  # R2: Pseudo_r2 function defined above
  r2_train <- cor(sitedata.train$logPM25, results.train)^2
  r2_val <- cor(sitedata.test$logPM25, results.val)^2
  
  r2_train_exp <- cor(sitedata.train$PM25, results.train_exp)^2
  r2_val_exp <- cor(sitedata.test$PM25, results.val)^2
  
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
  
  rmse_train_exp <- RMSE(sitedata.train$PM25, results.train_exp)
  rmse_val_exp <- RMSE(sitedata.test$PM25, results.val_exp)
  rmae_train_exp <- RMAE(sitedata.train$PM25, results.train_exp)
  rmae_val_exp <- RMAE(sitedata.test$PM25, results.val_exp)
  bias_train_exp <- bias(sitedata.train$PM25, results.train_exp)
  bias_val_exp <- bias(sitedata.test$PM25, results.val_exp)
  r2_train_exp <- r2(sitedata.train$PM25, results.train_exp)
  r2_val_exp <- r2(sitedata.test$PM25, results.val_exp)
  pmcc_train_exp <- PMCC(sitedata.train$PM25, results.train_exp)
  pmcc_val_exp <- PMCC(sitedata.test$PM25, results.val_exp)
  
  
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
                           "cov" = cov, # 71
                           "se_val_exp" = se_val_exp, # 11
                           "se_train_exp" = se_train_exp, # 16
                           "r2_val_exp" = r2_val_exp, # 21
                           "r2_train_exp" = r2_train_exp, # 26
                           "rmse_train_exp" = rmse_train_exp, # 31
                           "rmse_val_exp" = rmse_val_exp, # 36
                           "rmae_train_exp" = rmae_train_exp, # 41
                           "rmae_val_exp" = rmae_val_exp, # 46
                           "bias_train_exp" = bias_train_exp, # 51
                           "bias_val_exp" = bias_val_exp, # 56
                           "pmcc_train_exp" = pmcc_train_exp, # 61
                           "pmcc_val_exp" =  pmcc_val_exp) # 66
  
  
}    

results <- do.call(rbind, results_list)
save(results, file = "MODELS/MODEL/Outputs/M1a_exp_SpatialBlockCV.RData")

#load(file = "MODELS/Outputs/M1_New_PCM_NoInd_SpatialBlockCV.RData")

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




