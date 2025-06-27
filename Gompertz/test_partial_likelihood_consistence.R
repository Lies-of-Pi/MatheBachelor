# Test auf Konsistenz des Beta - Wertes durch Cox

library(parallel)
library(flexsurv)
library(ebmstate)
library(bindata)
library(Hmisc)

# Parameter um die richtigen files zu laden
n_state <- 9 # Länge der Kette
n_cov <- 5 # number of covariates ---> Mal ein ganz primitives Modell, damit ich garantiert genug daten habe
n_patients <- 10000
n_steps <- 100

datafile<-paste0("R_data/", n_patients, "_patients_", n_cov, "_covariates_",n_state,"states.RData")
model_file <-paste0("R_data/", n_patients, "_patients_", n_cov, "_covariates_",n_state,"states_true_model.RData")

# loading data
load(datafile)
load(model_file)

# get the true beta parameter
beta_true <- model_data$beta

# define beta function for the Cox Models
beta_cox <- function(trans){
  coxrfx_object$coefficients[seq(trans,length(coxrfx_object$coefficients),n_state - 1)]
}


# Get the data into right format
mstate.data.all <- data_long_format$mstate.data

#expand covariates
mstate.data.all<-mstate::expand.covs(mstate.data.all,covs =names(mstate.data.all)[-(1:8)])

MSE_matrix <- matrix(nrow = n_state - 1, ncol = n_steps)

for (i in 1:n_steps){
  print(i/n_steps)
  # get the right data
  data_mask <- mstate.data.all$id <= i * n_steps # Sum data_mask = Anzahl verwendeter Daten
  mstate.data <- mstate.data.all[data_mask, ]
  
  #argument 'Z' of coxrfx
  Z<-mstate.data[,-(1:(8+n_cov))] # Nur die Expanded covs
  Z$strata<-Z$trans<-mstate.data$trans # CoxRFX braucht noch die Daten über die Verwendeten Übergänge / Welche basline function welchen Übergang beschreibt
  
  #argument 'surv' of coxrfx
  surv<-survival::Surv(mstate.data$time,mstate.data$status)
  
  #argument 'groups' of coxrfx
  groups<-rep("unique_group", n_cov * (n_state - 1))
  
  #fit random effects model --> Dieses Modell setzt eine gemeinsame Verteilung der beta Werte voraus
  coxrfx_object<-CoxRFX(Z,surv,groups,max.iter = 600,tol = 0.0001,sigma.hat = "df")
  cat("coxrfx_concordance:",concordance(coxrfx_object)$concordance)
  
  # Schaue dir die Betragsquadrate der Vektoren an:
  for (trans in 1:(n_state - 1)){
    beta_true_trans <- beta_true[(1 + n_cov * (trans - 1)) : (n_cov * trans)]
    difference <- beta_true_trans - beta_cox(trans)
    MSE <- difference%*%difference
    MSE_matrix[trans, i] <- MSE
    print(paste0("Transition: ", trans, " MSE ", MSE))
  }
}

# Mache einen Plot um die Entwicklung der Differenzen zu Zeigen
n_data <- 1:n_steps
state_colors <- c("red", "green", "blue", "yellow", "pink", "cyan", "violet", "lightblue")
beta_true_length <- function(trans){
  beta_true_trans <- beta_true[(1 + n_cov * (trans - 1)) : (n_cov * trans)]
  out <- beta_true_trans %*% beta_true_trans
  out
}

pdf("Test_consistence_cox_5_cov.pdf")
plot(n_data,MSE_matrix[1,] / beta_true_length(1)[1] ,lwd=2,col=state_colors[1],lty=3, ylab="squared difference / squared length", xlab=paste0("used samples x ", n_steps), main="Test der Beta Konsistenz Cox", ylim = c(0, 1.5))
for (i in 2:(n_state - 1)){
  points(n_data,MSE_matrix[i,] / beta_true_length(i)[1],lwd=2,col=state_colors[i],lty=3)
}
legend(80, 1.5, legend=c(paste0("transition ", 1:(n_state - 1))), col=state_colors,lwd=2:2, lty=3:3, cex=0.8)
dev.off()