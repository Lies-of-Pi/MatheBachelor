# Teste die Varianz meiner Aussage
n <- 10000

# Generiere lambda Werte
lambda1_all_values <- c(0.2, 0.5, 1, 10, 25, 50, 100)
lambda2 <- seq(0.1, 100, 0.1)

var_values <- matrix(rep(0, length(lambda1_all_values) * length(lambda2)), ncol= length(lambda2) , nrow = length(lambda1_all_values))
for (i in 1:length(lambda1_all_values)){
  for (j in 1:length(lambda2)){
    # Bilde die Werte
    t1 <- rexp(n, lambda1_all_values[i])
    t2 <- rexp(n, lambda2[j])
    values <- log(t1/t2)
    var_values[i, j] = var(values)
  }
}
png("Empirical_proof_of_variance.png")
plot_colors <- c("red", "green", "blue", "yellow", "pink", "cyan", "violet", "black")
plot(lambda2, var_values[1,], ylim=c(0, 5), ylab="Var(log(t1/t2))", col=plot_colors[1], main="Empirical illustration of Variance statement")
for (i in 2:length(lambda1_all_values)){
  points(lambda2, var_values[i,], col=plot_colors[i])
}

# Mein theoretischer Wert
lines(lambda2, rep(3.269918, length(lambda2)), col=plot_colors[length(lambda1_all_values) + 1], lwd=5)
legend(0, 2, legend=c(c(paste0("lambda1 = ", lambda1_all_values), "Theoretical value")), col=plot_colors, lwd=4, lty=3:3, cex=0.8)
dev.off()