# Generierung der Bilder um die Daten zu prüfen
# Wichtig!!!!! Hier darf m nur Daten einer Covariate beinhalten

library(parallel)
library(flexsurv)
library(ebmstate)
library(bindata)
library(Hmisc)

# Parameter um die richtigen files zu laden
n_state <- 9 # Länge der Kette
n_cov <- 100 # number of covariates
n_patients <- 10000


# Parameter für die Plots / Empirische Auswertung
nTimepoints <- 10000

# Generierung der Filenamen
filenames <- paste0("StateProbTheory/state_",1:n_state, ".png")

datafile_m<-paste0("R_data/", n_patients, "_patients_", n_cov, "_covariates_",n_state,"states_m.RData")
model_file <-paste0("R_data/", n_patients, "_patients_", n_cov, "_covariates_",n_state,"states_true_model.RData")

# loading data
load(datafile_m)
load(model_file)

# Hier bekommen wir alle Daten der einen Covariate und die Covariatenvektor selbst
m <- m_data_for_test$m
covariate_vector <- m_data_for_test$covariate_vector

time_max <- 1.25 * max(m[, n_state-1])
time_vector<- seq(0,time_max,time_max / nTimepoints)


#A function that computes relative frequencies of each state at some time t
# Die Daten beschreiben jetzt die Übergangszeitpunkte, hier muss ich nur durchzählen
rel_freq<-function(state,t,m){
  # Erster Zustand:
  if (state == 1){
    sum(t < m[, state])/nrow(m)
  }
  else if (state != n_state){
    sum((t < m[, state]) & (t > m[, state-1]))/nrow(m) 
  }
  else if (state == n_state){
    sum((t > m[, n_state-1]))/nrow(m)
  }
}

# Hauptschleife: Berechne empirische und theoretische StateProb und Plotte diese
for (state in 1:n_state){
  print(paste0("state: ", state))
  
  # Bestimmung der empirischen Stateprob
  empirical_state_prob <- sapply(time_vector, function(t){rel_freq(state, t, m)})
  
  # Bestimmung meiner Stateprob
  theoretical_state_prob <- sapply(time_vector, function(t) {model_data$true_model(state, covariate_vector, model_data$a, model_data$b, model_data$beta, model_data$n_state, model_data$n_cov ,t)})
  
  # Plot der Ergebnisse
  png(filenames[state])
  plot(time_vector,theoretical_state_prob,lwd=2,col="green",lty=3, ylab=paste0("state ", state, " probability"), xlab="time", main=paste0("Gompertz with a = ", model_data$a, ", b = ", model_data$b[state]))
  lines(time_vector,empirical_state_prob,lwd=2,col="red",lty=3)
  legend(0, 0.1, legend=c("theoretical", "empirical"), col=c("green", "red"),lwd=2:2, lty=3:3, cex=0.8)
  dev.off()
}

