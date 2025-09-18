# Simulation Scenario 2
# Compare power of the MA and RB tests before and after covariate adjustment
# When adjusting, compare correct and incorrect model specification
rm(list = ls())
library(ggplot2)
library(lmtest)
library(sandwich)

# Number of simulations, set high for final run before submission
nsims <- 1000

# Initialize result arrays
pate_true <- cate_true <- array(NA, dim = c(nsims, 8, 1))
diff_pval1 <- diff_pval2 <- diff_pval3 <- diff_pval4 <- 
  diff_pval5 <- diff_pval6 <- array(NA, dim = c(nsims, 8, 1))

# Set random seed
set.seed(1234)

# Simulation
# Loop over 3 number of clusters
for (j in seq(30, 100, by = 10)) {
  
  # Loop over null effect size
  for (k in 5:5) {
    
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
      
      # Informative cluster level covariate
      CC <- rnorm(num_clusters, 0, 1)
      CC_full <- rep(CC, cluster_sizes)
      
      # Simulate potential outcomes with ICC = 0.1
      ICC <- 0.1
      Y_resvar <- 1
      alpha0_var <- Y_resvar * ICC / (1 - ICC)
      alpha0 <- rnorm(num_clusters, 0, sqrt(alpha0_var))
      Y0 <- rep(NA, n)
      Y0[sizecheck == 0] <- rnorm(sum(1-sizecheck), 0 + 1*CC_full,
                                  Y_resvar) +
        rep(alpha0[cluster_sizes <= 50], cluster_sizes[cluster_sizes <= 50])
      Y0[sizecheck == 1] <- rnorm(sum(sizecheck), 1 + 1*CC_full,
                                  Y_resvar) +
        rep(alpha0[cluster_sizes > 50], cluster_sizes[cluster_sizes > 50])
      
      Y1 <- rep(NA, n)
      Y1[sizecheck == 0] <- rnorm(sum(1-sizecheck), 0.5 + 1*CC_full, Y_resvar) +
        rep(alpha0[cluster_sizes <= 50], cluster_sizes[cluster_sizes <= 50])
      Y1[sizecheck == 1] <- rnorm(sum(sizecheck),
                                  1.5 + 1*(k-5) + 1*CC_full, Y_resvar) +
        rep(alpha0[cluster_sizes > 50], cluster_sizes[cluster_sizes > 50])
      
      # Exact 1:1 treatment allocation
      trtprob <- 0.5
      trt <- sample(1:num_clusters, num_clusters*trtprob)
      Z <- rep(0, n)
      Z[clusternum %in% trt] <- 1
      
      # Observed outcome and true values
      Y <- Z*Y1 + (1-Z)*Y0
      pate_true[i, floor(j/10 - 2), k-4] <- mean(Y1 - Y0)
      cate_true[i, floor(j/10 - 2), k-4] <- sum((Y1 - Y0) / cluster_sizes_full) /
        num_clusters
      
      # Proposed model assisted tests
      Yvec <- tapply(Y, clusternum, mean)
      Zvec <- tapply(Z, clusternum, mean)
      w <- cluster_sizes / n - 1 / num_clusters
      wYvec <- w*Yvec*num_clusters
      cs_cent <- 1*(cluster_sizes > 50) - mean(1*(cluster_sizes > 50))
      
      # Unadjusted MAT
      diff_pval1[i, floor(j/10 - 2), k-4] <-
        coeftest(lm(wYvec ~ Zvec), vcov = vcovHC(lm(wYvec ~ Zvec),
                                                 type = "HC0"))[2, 4]
      
      # Adjusted MAT for cluster size and cluster covariate
      diff_pval2[i, floor(j/10 - 2), k-4] <-
        coeftest(lm(wYvec ~ Zvec*CC + Zvec*cs_cent),
                 vcov = vcovHC(lm(wYvec ~ Zvec*CC + Zvec*cs_cent),
                               type = "HC0"))[2, 4]
      
      # Misspecified adjusted MAT
      CC2 <- CC^2
      diff_pval3[i, floor(j/10 - 2), k-4] <-
        coeftest(lm(wYvec ~ Zvec*cs_cent + Zvec*CC2),
                 vcov = vcovHC(lm(wYvec ~ Zvec*cs_cent + Zvec*CC2),
                               type = "HC0"))[2, 4]
      
      # Randomization based tests done similarly
      numperms <- 200
      statvec1 <- statvec2 <- statvec3 <- rep(NA, numperms)
      for (perm in 1:numperms) {
        
        trtperm <- sample(1:num_clusters, num_clusters*trtprob)
        Zperm <- rep(0, n)
        Zperm[clusternum %in% trtperm] <- 1
        Zvecperm <- tapply(Zperm, clusternum, mean)
        
        statvec1[perm] <- coeftest(lm(wYvec ~ Zvecperm),
                                   vcov = vcovHC(lm(wYvec ~ Zvecperm),
                                                 type = "HC0"))[2, 3]
        statvec2[perm] <-
          coeftest(lm(wYvec ~ Zvecperm*CC + Zvecperm*cs_cent),
                   vcov = vcovHC(lm(wYvec ~ Zvecperm*CC + Zvecperm*cs_cent),
                                 type = "HC0"))[2, 3]
        statvec3[perm] <-
          coeftest(lm(wYvec ~ Zvecperm*cs_cent + Zvecperm*CC2),
                   vcov = vcovHC(lm(wYvec ~ Zvecperm*cs_cent + Zvecperm*CC2),
                                 type = "HC0"))[2, 3]
        
      }
      
      diff_pval4[i, floor(j/10 - 2), k-4] <-
        mean(abs(statvec1) >= abs(coeftest(lm(wYvec ~ Zvec),
                                           vcov = vcovHC(lm(wYvec ~ Zvec),
                                                         type = "HC0"))[2, 3]))
      diff_pval5[i, floor(j/10 - 2), k-4] <-
        mean(abs(statvec2) >=
               abs(coeftest(lm(wYvec ~ Zvec*CC + Zvec*cs_cent),
                            vcov = vcovHC(lm(wYvec ~ Zvec*CC + Zvec*cs_cent),
                                          type = "HC0"))[2, 3]))
      diff_pval6[i, floor(j/10 - 2), k-4] <-
        mean(abs(statvec3) >=
               abs(coeftest(lm(wYvec ~ Zvec*cs_cent + Zvec*CC2),
                            vcov = vcovHC(lm(wYvec ~ Zvec*cs_cent + Zvec*CC2),
                                          type = "HC0"))[2, 3]))
      
    }
    
  }
  
}


# Make results figure
plotdat <- data.frame("Method_Type" = rep(c("Model-assisted Test",
                                            "Randomization-based test"), each = 24),
                      "Method" = c(rep("MAT unadjusted", 8),
                                   rep("MAT adjusted", 8),
                                   rep("MAT misspecified", 8),
                                   rep("RBT unadjusted", 8),
                                   rep("RBT adjusted", 8),
                                   rep("RBT misspecified", 8)),
                      "Num_Clusters" = rep(rep(seq(30, 100, by = 10), each = 1), 6),
                      "T1E" = c(apply(diff_pval1[, , ], 2, function(x) mean(x < 0.05)),
                                apply(diff_pval2[, , ], 2, function(x) mean(x < 0.05)),
                                apply(diff_pval3[, , ], 2, function(x) mean(x < 0.05)),
                                apply(diff_pval4[, , ], 2, function(x) mean(x < 0.05)),
                                apply(diff_pval5[, , ], 2, function(x) mean(x < 0.05)),
                                apply(diff_pval6[, , ], 2, function(x) mean(x < 0.05))))

# Make figure
ggplot(dat = plotdat, aes(x = Num_Clusters, y = T1E,
                          color = Method, shape = Method, group = Method)) +
  geom_line() +
  geom_point() +
  #facet_grid(~Method) +
  geom_hline(yintercept = 0.05) +
  geom_hline(yintercept = 0.05 + 1.96*sqrt(0.05*0.95/nsims), lty = 2) +
  geom_hline(yintercept = 0.05 - 1.96*sqrt(0.05*0.95/nsims), lty = 2) +
  xlab("Number of Clusters") +
  ylab("Type I Error") +
  theme_bw()

# Save results
save(pate_true, cate_true, diff_pval1, diff_pval2, diff_pval3,
     diff_pval4, diff_pval5, diff_pval6,
     file = "~/OneDrive_VUMC/Projects/ICS/Output/sims-scenario2-t1e-appendix.RData")

