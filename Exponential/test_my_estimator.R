# Test des Konsistenzverhaltens meines Schätzers

library(parallel)
library(flexsurv)
library(ebmstate)
library(bindata)
library(Hmisc)

# Parameter um die richtigen files zu laden
n_state <- 7 # Länge der Kette
n_cov <- 5 # number of covariates ---> Mal ein ganz primitives Modell, damit ich garantiert genug Daten habe
n_patients <- 100000
n_steps <- 100

datafile<-paste0("R_data/", n_patients, "_patients_", n_cov, "_covariates_",n_state,"states_exp_m.RData")
model_file <-paste0("R_data/", n_patients, "_patients_", n_cov, "_covariates_",n_state,"states_true_model_exp.RData")

# loading data
load(datafile)
load(model_file)

# get the true beta parameter
beta_true <- model_data$beta

# Definiere meinen Schätzer
beta_estimator <- function(data){
  beta_hat <- c()
  n <- length(data[, 1])
  
  # 1. Berechne die Varianzen der Covariaten und bilde die Vorfaktormatrix:
  A <- diag(n_cov)
  for (i in 1:n_cov){
    cov_var <- var(data[, n_state - 1 + i])
    if(cov_var == 0) stop ("Varianz einer Covariate ist 0!!!!")
    A[i, i] <- 1/cov_var
  }
  
  X <- data[, n_state:(n_state + n_cov - 1)]
  times <- data[, 1:(n_state - 1)]
  
  #2. Schleife über alle Übergänge
  for (trans in 1:(n_state - 1)){
    t <- times[, trans]
    
    # 3. Berechne den Vektor über (mean(ln(t)) * mean(x) - mean(ln(t) * x))
    mean_ln_t <- mean(log(t))
    mean_x <- colMeans(X)
    mean_ln_t_x <- colMeans(log(t) * X)
    
    vector <- (mean_ln_t * mean_x - mean_ln_t_x)
    
    # 4. Berechne beta:
    beta_trans <- A %*% vector
    
    # 5. Speichere den Vektor
    beta_hat <- c(beta_hat, beta_trans)
    
  }
  beta_hat
}



n_points <- n_patients / 100
MSE_matrix <- matrix(nrow = n_state - 1, ncol = n_points)

for (i in 1:n_points){
  print(i/n_points)
  # get the right data
  data <- m[1:(n_steps*i), ]
  
  # Schätze beta
  beta_hat <- beta_estimator(data)
  
  # Schaue dir die Betragsquadrate der Vektoren an:
  for (trans in 1:(n_state - 1)){
    beta_true_trans <- beta_true[(1 + n_cov * (trans - 1)) : (n_cov * trans)]
    beta_hat_trans <- beta_hat[(1 + n_cov * (trans - 1)) : (n_cov * trans)]
    difference <- beta_true_trans - beta_hat_trans
    MSE <- difference%*%difference
    MSE_matrix[trans, i] <- MSE
    print(paste0("Transition: ", trans, " MSE ", MSE))
  }
}

# Mache einen Plot um die Entwicklung der Differenzen zu Zeigen
n_data <- 1:n_points
state_colors <- c("red", "green", "blue", "yellow", "pink", "cyan", "violet", "lightblue")
beta_true_length <- function(trans){
  beta_true_trans <- beta_true[(1 + n_cov * (trans - 1)) : (n_cov * trans)]
  out <- beta_true_trans %*% beta_true_trans
  out
}

x <-log(100 * n_data)
y <-log(MSE_matrix[1,] / beta_true_length(1)[1])
pdf("Test_consistence_exp_estimator.pdf")
plot(x,y ,lwd=2,col=state_colors[1],lty=3, ylab="log(squared difference / squared length)", xlab="log(used samples)", main="Test my estimator")
grid()
for (i in 2:(n_state - 1)){
  points(x, log(MSE_matrix[i,] / beta_true_length(i)[1]),lwd=2,col=state_colors[i],lty=3)
}
legend(10, -2, legend=c(paste0("transition ", 1:(n_state - 1))), col=state_colors[1:(n_state - 1)],lwd=2:2, lty=3:3, cex=0.8)
dev.off()

png("Test_consistence_exp_estimator.png")
plot(x,y ,lwd=2,col=state_colors[1],lty=3, ylab="log(squared difference / squared length)", xlab="log(used samples)", main="Test my estimator")
for (i in 2:(n_state - 1)){
  points(x,log(MSE_matrix[i,] / beta_true_length(i)[1]),lwd=2,col=state_colors[i],lty=3)
}
legend(10, -2, legend=c(paste0("transition ", 1:(n_state - 1))), col=state_colors[1:(n_state - 1)],lwd=2:2, lty=3:3, cex=0.8)
dev.off()