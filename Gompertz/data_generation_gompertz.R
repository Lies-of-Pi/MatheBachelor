# Data generation for linear chain with my model
# Dieser Code generiert die nötigen Daten und speichert sie in einer csv Datei ab.
# Die Daten des Wahren modells sollen in einer Rdata Datei abgespeichert werden.
# Der Code Speichert außerdem einen Histogramm file um die Streuung der Daten mit den patienten zu veranschaulichen

library(parallel)
library(flexsurv)
library(ebmstate)
library(bindata)
library(Hmisc)

# Step 1 Define Parameters ###############################################################
# Definiere Parameter der Datengenerierung 
n_state <- 9 # Länge der Kette
n_cov <- 100 # number of covariates
n_patients <- 10000 # number of patients
censoring_level <- 0.3 # we want a maximum of censoring_level of censored data
just_one_covariate <- FALSE # generate data for n_patients time for the same patient

a <- 1 #0.1 # Parameter fuer die Gompertz Verteilung
epsilon <- 0.001 
b_1 <- exp(runif(1, min=-3, -1))
lambda_censoring <- 0.0004

datafile<-paste0("R_data/", n_patients, "_patients_", n_cov, "_covariates_",n_state,"states.RData")
datafilecsv <- paste0("R_data/", n_patients, "_patients_", n_cov, "_covariates_",n_state,"states.csv")
datafile_m<-paste0("R_data/", n_patients, "_patients_", n_cov, "_covariates_",n_state,"states_m.RData")
model_file <-paste0("R_data/", n_patients, "_patients_", n_cov, "_covariates_",n_state,"states_true_model.RData")
hist_file <- paste0("R_data/", n_patients, "_patients_", n_cov, "_covariates_",n_state,"states_histogram.png")

# Define automatic parameters for the generation
covariate_names<-paste0("Cov",1:n_cov)
nTrans <- n_state - 1
filenames <- paste0("StateProbTheory/state_",1:n_state, ".png")
covariate_names<-paste0("Cov",1:n_cov)
tmat_transitions <- list(c())
for (i in n_state:2){
  tmat_transitions <- append(tmat_transitions, c(i))
}
tmat<-mstate::transMat(x=rev(tmat_transitions),names=c(paste0("state", 1:n_state)))

# Step 2 / 3 define covariates ####################################################
marg_probs<-runif(n = length(covariate_names),min = 0.05,max = 0.3)

if (just_one_covariate){
  covariate_vector <- rmvbin(1,margprob = marg_probs)
  covariate_matrix <- matrix(covariate_vector, ncol = n_cov)
  for (i in 2:n_patients){
    covariate_matrix <- rbind(covariate_matrix, covariate_vector)
  }
}
if (!just_one_covariate){
  covariate_matrix<-rmvbin(n_patients,margprob = marg_probs) #cov matrix 
}
colnames(covariate_matrix)<-covariate_names

# Step 4 Erzeugung der Modelparameter des wahren Modells ##########################################
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
  
  out <- min(min_set) * exp(- y_i_dach %*% beta_i) * epsilon / (1 - epsilon)
  out
}

# Generate b:
b <- array(b_1)
for (i in 2:nTrans){
  b <- append(b, b_generation(i, b, epsilon))
}

# Step 5: Erzeugung der Daten #####################################################
risk_trans <- function(i){
  out <- exp(covariate_matrix%*%beta[(1+length(covariate_names)*(i - 1)):(length(covariate_names)*i)])
  out
}
rates <- sapply(1:(n_state-1), function(j){risk_trans(j) * b[j]})

# Generate Data
m <- matrix(sapply(1:(n_state - 1), function(i){flexsurv::rgompertz(n_patients, shape=a, rate = rates[, i]) }), ncol = n_state - 1)
m<-cbind(m,rexp(n_patients,lambda_censoring)) 
m<-cbind(m,(m[,n_state]<m[,n_state - 1]))
colnames(m)<-c(paste0("T_", 1:(n_state - 1)), "cens_time", "cens=1")
m<-as.data.frame(m)
View(m)

# falls wir nur eine covariate generiert haben, vergleich mit dem Theoretischen Modell
if (just_one_covariate){
  m_data_for_test <- list(m=m, covariate_vector=covariate_matrix[1,])
  save(m_data_for_test, file=datafile_m)
}

# Sehr wichtiges Maß, Anteil censoring:
censoring_ratio <- sum(m[, n_state + 1])/n_patients
print("Censoring ratio:")
print(censoring_ratio)
if(censoring_ratio > censoring_level) stop('Too much censoring, choose smaller censoring lambda')


# Step 6: Transform Data into long format ##############################
#convert the data to long format
mstate.data<-data.frame() # Daten werden in anderes format gebracht

for(i in 1:nrow(m)){ 
  # Generate transition history for each patient:
  status <- 1
  for (j in 1:(n_state - 1)){
    # check if transition failed
    if (status == 0) break
    id <- i
    from <- j # sind gerade in Zustand j
    to <- j + 1
    trans <- j
    if (j == 1) Tstart <- 0
    else Tstart <- m[i, j-1] # Wann sind wir nach j gekommen?
    Tstop <- min(m[i, j],m$cens_time[i]) # findet der Übergang vor der Censierung statt?
    time <- Tstop - Tstart
    status <- as.numeric(m$cens_time[i]>m[i, j])
    mstate.data<-rbind(mstate.data,data.frame(id=id,from=from,to=to,trans=trans,Tstart=Tstart,Tstop=Tstop,time=time,status=status))
    
  }
}
#add covariates
mstate.data<-cbind(mstate.data,covariate_matrix[mstate.data$id,])
View(mstate.data)

# Step 7: Export Data ####################################################################
#attributes and class
class(mstate.data)<-c("data.frame","msdata")
attr(mstate.data,"trans")<-tmat

# Speicherung der Daten
data_long_format <- list(mstate.data=mstate.data, beta=beta, a=a, rate=rates, b=b)
save(data_long_format, file=datafile)

# Speichere die Matrix besser anders ab
write.csv(mstate.data, datafilecsv, row.names = TRUE)


# Step 8: Creating true theoretical model ##############################################
true_model <- function(state, covariate_vector, a, b, beta, n_state, n_cov, t){
  # Liefert die Besetzungswkt zur Covariate cov zum Zeitpunkt t, macht nur sinn
  # wenn covariate_vector ein vector ist.
  
  # Definition der Konstanten
  # c-Vektor mit allen c_i elementen
  c <- sapply(1:(n_state - 1), function(i){b[i]/a * exp(covariate_vector%*%beta[(1+n_cov*(i-1)):(n_cov*i)])})
  
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
    out <- sum(sapply((1:state), function(j){ D[j] * exp(- c[j] * exp(a*t)) * Product_matrix[state, j]}))
  }
  else{
    D_n <- sum(sapply(1:(n_state-1), function(j){ D[j] * exp(- c[j]) * Product_matrix[n_state-1, j] *(c[n_state - 1]/c[j])}))
    out <- D_n - sum(sapply((1:(n_state - 1)), function(j){ D[j] * exp(- c[j] * exp(a*t)) * Product_matrix[n_state - 1, j]  *(c[n_state - 1]/c[j])}))
  }
  out
}

# save the true model
model_data <- list(true_model=true_model, beta=beta, a=a, b=b, n_state=n_state, n_cov=n_cov)
save(model_data, file=model_file)


# Step 9: Creating histogram plot of the Time Data
#hist_colors <- sample(colors())
png(hist_file)
hist_colors <- c("red", "blue", "green", "yellow", "pink", "cyan", "magenta", "lightblue", "lightgreen", "orange")
xmin <- 0
xmax <- 1.1 * max(m[,n_state - 1])
hist(m[, 1], xlab="time", col=hist_colors[1], xlim=c(xmin, xmax), freq=FALSE, main="Data distribution")
for (state in 2:(n_state - 1)){
  hist(m[, state], xlab="time", col=hist_colors[state], add=TRUE, xlim=c(xmin, xmax), freq=FALSE)
}

legend("topright", colnames(m)[1:(n_state - 1)], 
       bty = "n", fill=hist_colors[1:(n_state - 1)])
dev.off()


