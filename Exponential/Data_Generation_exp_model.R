#Data generation for Exp distribution

library(parallel)
library(flexsurv)
library(ebmstate)
library(bindata)
library(Hmisc)

# Step 1 Define Parameters ###############################################################
# Definiere Parameter der Datengenerierung 
n_state <- 7 # Länge der Kette
n_cov <- 5 # number of covariates
n_patients <- 100000 # number of patients

epsilon <- 0.001 
b_factor <- 10**2 ########################## Ändern falls es schief geht
b_1 <- exp(runif(1, min=-3, -1))

datafile_m<-paste0("R_data/", n_patients, "_patients_", n_cov, "_covariates_",n_state,"states_exp_m.RData")
model_file <-paste0("R_data/", n_patients, "_patients_", n_cov, "_covariates_",n_state,"states_true_model_exp.RData")

# Define automatic parameters for the generation
covariate_names<-paste0("Cov",1:n_cov)
nTrans <- n_state - 1
filenames <- paste0("StateProbTheory/state_",1:n_state, ".png")
covariate_names<-paste0("Cov",1:n_cov)


# Define covariates ####################################################
marg_probs<-runif(n = length(covariate_names),min = 0.05,max = 0.3)

covariate_matrix<-rmvbin(n_patients,margprob = marg_probs) #cov matrix 
colnames(covariate_matrix)<-covariate_names

# Erzeugung der Modelparameter des wahren Modells ##########################################
# Beta Parameter:
nParam<-nTrans*length(covariate_names)
param_fun<-function(nParam){
  out<-sqrt(10/length(covariate_names))*rnorm(n = nParam,mean = 0,sd = 0.65)
  out 
}
beta<-param_fun(nParam) # Beta parameter

# b parameter:
b_generation <- function(i, b, epsilon){
  # Wir wollen b für i berechnen, i > 1
  
  # Berechne den min Faktor
  min_set <- c()
  for (k in 1:(i - 1)){
    beta_k <- beta[(1+length(covariate_names)*(k - 1)):(length(covariate_names)*k)]
    y_k <- (beta_k < 0)
    min_set <- c(min_set, b[k] * exp(y_k %*% beta_k))
  }
  beta_i <- beta[(1+length(covariate_names)*(i - 1)):(length(covariate_names)*i)]
  y_i_dach <- (beta_i > 0)
  
  out <- min(min_set) * exp(- y_i_dach %*% beta_i) * epsilon / (1 - epsilon) * b_factor
  out
}

# Generate b:
b <- array(b_1)
for (i in 2:nTrans){
  b <- append(b, b_generation(i, b, epsilon))
}

# Erzeugung der Daten #####################################################
risk_trans <- function(i){
  out <- exp(covariate_matrix%*%beta[(1+length(covariate_names)*(i - 1)):(length(covariate_names)*i)])
  out
}
rates <- sapply(1:(n_state-1), function(j){risk_trans(j) * b[j]})

# Generate Data
m <- matrix(sapply(1:(n_state - 1), function(i){rexp(n_patients, rate = rates[, i]) }), ncol = n_state - 1)
colnames(m)<-c(paste0("T_", 1:(n_state - 1)))
m<-as.data.frame(m)

# add covariates
m<-cbind(m,covariate_matrix)
View(m)

#  Export Data ####################################################################
save(m, file=datafile_m)



# Creating true theoretical model ##############################################
true_model <- function(state, covariate_vector, b, beta, n_state, n_cov, t){
  # Liefert die Besetzungswkt zur Covariate cov zum Zeitpunkt t, macht nur sinn
  # wenn covariate_vector ein vector ist.
  
  # Definition der Konstanten
  # c-Vektor mit allen c_i elementen
  c <- sapply(1:(n_state - 1), function(i){b[i] * exp(covariate_vector%*%beta[(1+n_cov*(i-1)):(n_cov*i)])})
  
  # Produkt Faktor als Matrix, an der stelle i, j soll der Wert des Produkts stehen für j+1 bis i stehen, für i = 1 ist das Produkt 1
  Product_matrix <- matrix(1, nrow = n_state - 1, ncol = n_state - 1)
  for (i in 2:(n_state - 1)){
    for (j in 1:(i - 1)){
      Product_matrix[i, j] <- prod(c[j:(i - 1)]) / (prod(c[(j + 1):i] - c[j]))
    }
  }
  
  # Vektor mit den Werten für D:
  D <- matrix(exp(c[1]), nrow = 1, ncol = 1)
  for (i in 2:(n_state - 1)){
    D[i] <- - sum(sapply(1:(i-1), function(j){ D[j] * exp(c[i] - c[j]) * Product_matrix[i, j]}))
  }
  
  # Eigentlicher Funktionsaufruf
  if (state < n_state){
    out <- sum(sapply((1:state), function(j){ D[j] * exp(- c[j]) * Product_matrix[state, j]}))
  }
  else{
    D_n <- sum(sapply(1:(n_state-1), function(j){ D[j] * exp(- c[j]) * Product_matrix[n_state-1, j] *(c[n_state - 1]/c[j])}))
    out <- D_n - sum(sapply((1:(n_state - 1)), function(j){ D[j] * exp(- c[j]) * Product_matrix[n_state - 1, j]  *(c[n_state - 1]/c[j])}))
  }
  out
}

# save the true model
model_data <- list(true_model=true_model, beta=beta, b=b, n_state=n_state, n_cov=n_cov)
save(model_data, file=model_file)


