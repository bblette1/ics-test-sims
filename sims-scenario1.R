# Simulation Scenario 1
# Compare proposed unadjusted methods to naive method
rm(list = ls())
library(ggplot2)
library(geepack)
library(lmtest)
library(sandwich)

# Number of simulations, set high for final run before submission
nsims <- 2000

# Initialize result arrays
pate_true <- cate_true <- array(NA, dim = c(nsims, 3, 9))
diff_pval1 <- diff_pval2 <- diff_pval3 <-
  diff_pval4 <- diff_pval5 <- array(NA, dim = c(nsims, 3, 9))

# Set random seed
set.seed(1234)

# Simulation
# Loop over 3 number of clusters
for (j in c(30, 60, 100)) {
  
  # Loop over 9 effect sizes
  for (k in 1:9) {
    
    # Simulate nsims iterations
    for (i in 1:nsims) {
      
      # Cluster sizes formerly Poisson with mean 50, adjusted to get more spread
      # without implausible sizes
      num_clusters <- j
      cluster_sizes <- round(c(runif(num_clusters/2, 20, 50),
                               runif(num_clusters/2, 50, 80)))
      n <- sum(cluster_sizes)
      cluster_sizes_full <- rep(cluster_sizes, cluster_sizes)
      clusternum <- rep(1:num_clusters, cluster_sizes)
      
      # Clusters larger than 50 will have a different treatment effect (across k)
      sizecheck <- 1*(cluster_sizes_full > 50)
      
      # Simulate potential outcomes with ICC = 0.1
      ICC <- 0.1
      Y_resvar <- 1
      alpha0_var <- Y_resvar * ICC / (1 - ICC)
      alpha0 <- rnorm(num_clusters, 0, sqrt(alpha0_var))
      Y0 <- rep(NA, n)
      Y0[sizecheck == 0] <- rnorm(sum(1-sizecheck), 0, Y_resvar) +
        rep(alpha0[cluster_sizes <= 50], cluster_sizes[cluster_sizes <= 50])
      Y0[sizecheck == 1] <- rnorm(sum(sizecheck), 1, Y_resvar) +
        rep(alpha0[cluster_sizes > 50], cluster_sizes[cluster_sizes > 50])
      
      Y1 <- rep(NA, n)
      Y1[sizecheck == 0] <- rnorm(sum(1-sizecheck), 0.5, Y_resvar) +
        rep(alpha0[cluster_sizes <= 50], cluster_sizes[cluster_sizes <= 50])
      Y1[sizecheck == 1] <- rnorm(sum(sizecheck), 1.5 + 1*(k-5), Y_resvar) +
        rep(alpha0[cluster_sizes > 50], cluster_sizes[cluster_sizes > 50])
      
      # Exact 1:1 treatment allocation
      trtprob <- 0.5
      trt <- sample(1:num_clusters, num_clusters*trtprob)
      Z <- rep(0, n)
      Z[clusternum %in% trt] <- 1
      
      # Observed outcome and check true differences (or switch to DGM with known)
      Y <- Z*Y1 + (1-Z)*Y0
      pate_true[i, floor(j/30), k] <- mean(Y1 - Y0)
      cate_true[i, floor(j/30), k] <- sum((Y1 - Y0) / cluster_sizes_full) / num_clusters
      
      # Model-based tests (naive, correct specified, misspecified)
      cluster_sizes_full_50 <- 1*(cluster_sizes_full > 50)
      diff_pval1[i, floor(j/30), k] <-
        summary(geeglm(Y ~ Z + cluster_sizes_full_50, corstr = "independence",
             id = clusternum))$coefficients[3, 4]
      
      diff_pval2[i, floor(j/30), k] <-
        summary(geeglm(Y ~ Z*cluster_sizes_full_50, corstr = "independence",
                       id = clusternum))$coefficients[4, 4]
      
      cluster_sizes_full_log <- log(cluster_sizes_full)
      diff_pval3[i, floor(j/30), k] <-
        summary(geeglm(Y ~ Z*cluster_sizes_full_log, corstr = "independence",
                       id = clusternum))$coefficients[4, 4]
      
      # Proposed model assisted test
      Yvec <- tapply(Y, clusternum, mean)
      Zvec <- tapply(Z, clusternum, mean)
      w <- cluster_sizes / n - 1 / num_clusters
      wYvec <- w*Yvec*num_clusters
      diff_pval4[i, floor(j/30), k] <- coeftest(lm(wYvec ~ Zvec),
                                                 vcov = vcovHC(lm(wYvec ~ Zvec),
                                                               type = "HC0"))[2, 4]
      
      # Proposed randomization based test
      # Numperms should be large enough to preserve accuracy without too much
      # computational burden. Make high for final run before submission
      numperms <- 200
      statvec <- rep(NA, numperms)
      for (perm in 1:numperms) {
        
        trtperm <- sample(1:num_clusters, num_clusters*trtprob)
        Zperm <- rep(0, n)
        Zperm[clusternum %in% trtperm] <- 1
        Zvecperm <- tapply(Zperm, clusternum, mean)
        statvec[perm] <- coeftest(lm(wYvec ~ Zvecperm),
                                  vcov = vcovHC(lm(wYvec ~ Zvecperm),
                                                type = "HC0"))[2, 3]
        
      }

      diff_pval5[i, floor(j/30), k] <-
        mean(abs(statvec) >= abs(coeftest(lm(wYvec ~ Zvec),
                                          vcov = vcovHC(lm(wYvec ~ Zvec),
                                                        type = "HC0"))[2, 3]))
      
    }
    
  }
  
}


# Make results figure
plotdat <- data.frame("Method" = c(rep("Model-Based Naive", 27),
                                   rep("Model-Based", 27),
                                   rep("Model-Based Misspecified", 27),
                                   rep("Model-Assisted", 27),
                                   rep("Randomization-Based", 27)),
                      "Num_Clusters" = rep(rep(c(" 30 Clusters", " 60 Clusters",
                                                 "100 Clusters"), each = 9), 5),
                      "Effect_Size" = rep(rep(1:9, 3), 5),
                      "Power" = c(colMeans(diff_pval1[, 1, ] < 0.05),
                                  colMeans(diff_pval1[, 2, ] < 0.05),
                                  colMeans(diff_pval1[, 3, ] < 0.05),
                                  colMeans(diff_pval2[, 1, ] < 0.05),
                                  colMeans(diff_pval2[, 2, ] < 0.05),
                                  colMeans(diff_pval2[, 3, ] < 0.05),
                                  colMeans(diff_pval3[, 1, ] < 0.05),
                                  colMeans(diff_pval3[, 2, ] < 0.05),
                                  colMeans(diff_pval3[, 3, ] < 0.05),
                                  colMeans(diff_pval4[, 1, ] < 0.05),
                                  colMeans(diff_pval4[, 2, ] < 0.05),
                                  colMeans(diff_pval4[, 3, ] < 0.05),
                                  colMeans(diff_pval5[, 1, ] < 0.05),
                                  colMeans(diff_pval5[, 2, ] < 0.05),
                                  colMeans(diff_pval5[, 3, ] < 0.05)))

# Plot with approximate x-axis, mark in footnote or switch to DGM with known
ggplot(dat = plotdat, aes(x = (Effect_Size-5)/6.667, y = Power,
                               shape = Method, color = Method)) +
  geom_line() +
  geom_point() +
  facet_wrap(~Num_Clusters) +
  geom_hline(yintercept = 0.05) +
  geom_hline(yintercept = 0.05 + 1.96*sqrt(0.05*0.95/nsims), lty = 2) +
  geom_hline(yintercept = 0.05 - 1.96*sqrt(0.05*0.95/nsims), lty = 2) +
  xlab("i-ATE \u2013 c-ATE") +
  ylab("Power / Type I Error") +
  theme_bw()

# Save results
save(pate_true, cate_true, diff_pval1, diff_pval2, diff_pval3,
     diff_pval4, diff_pval5,
     file = "~/OneDrive_VUMC/Projects/ICS/Output/sims-scenario1.RData")
