# Simulation Scenario 3 (formerly scenario 4)
# Investigate the impact of an intermediate ICS test on primary bias/efficiency
# Assume that two analysts want to estimate the PATE
# Analyst 1 uses GEE with independence correlation matrix
# Analyst 2 uses an ICS test and then chooses an estimator based on the result
# Vary degree of ICS
rm(list = ls())
library(geepack)
library(ggplot2)
library(lmtest)
library(lme4)
library(sandwich)
library(ggpubr)
library(matrixStats)

# Number of simulations, set high for final run before submission
nsims <- 2000

# Initialize result arrays
pate_true <- cate_true <- array(NA, dim = c(nsims, 3, 9))
pate_bias1 <- pate_bias2 <- pate_bias3 <- array(NA, dim = c(nsims, 3, 9))
pate_se1 <- pate_se2 <- pate_se3 <- array(NA, dim = c(nsims, 3, 9))
pate_pval1 <- pate_pval2 <- pate_pval3 <- array(NA, dim = c(nsims, 3, 9))

# Set random seed
set.seed(1234)

# Simulation
# Loop over 1 number of clusters
for (j in c(100)) {
  
  # Loop over 9 effect sizes
  for (k in 1:9) {
    
    # Simulate nsims iterations
    for (i in 1:nsims) {
      
      # Cluster sizes formerly Poisson with mean 50, adjusted to get more spread
      # without implausible sizes
      # With mean 50 cluster size, analyst 3 would almost always become analyst 1
      # More different results with mean 20 cluster size
      num_clusters <- j
      #cluster_sizes <- round(c(runif(num_clusters/2, 20, 50),
                               #runif(num_clusters/2, 50, 80)))
      cluster_sizes <- round(c(runif(num_clusters/2, 5, 20),
                               runif(num_clusters/2, 20, 35)))
      n <- sum(cluster_sizes)
      cluster_sizes_full <- rep(cluster_sizes, cluster_sizes)
      clusternum <- rep(1:num_clusters, cluster_sizes)
      
      # Clusters larger than 50 will have a different treatment effect (across k)
      #sizecheck <- 1*(cluster_sizes_full > 50)
      sizecheck <- 1*(cluster_sizes_full > 20)
      
      # Cluster level covariate
      CC <- rnorm(num_clusters, 0, 1)
      CC_full <- rep(CC, cluster_sizes)
      
      # Simulate potential outcomes with ICC = 0.1
      ICC <- 0.1
      Y_resvar <- 1
      alpha0_var <- Y_resvar * ICC / (1 - ICC)
      alpha0 <- rnorm(num_clusters, 0, sqrt(alpha0_var))
      Y0 <- rep(NA, n)
      #Y0[sizecheck == 0] <- rnorm(sum(1-sizecheck), 0, Y_resvar) +
        #rep(alpha0[cluster_sizes <= 50], cluster_sizes[cluster_sizes <= 50])
      #Y0[sizecheck == 1] <- rnorm(sum(sizecheck), 1, Y_resvar) +
        #rep(alpha0[cluster_sizes > 50], cluster_sizes[cluster_sizes > 50])
      Y0[sizecheck == 0] <- rnorm(sum(1-sizecheck), 0 + 1*CC_full, Y_resvar) +
        rep(alpha0[cluster_sizes <= 20], cluster_sizes[cluster_sizes <= 20])
      Y0[sizecheck == 1] <- rnorm(sum(sizecheck), 1 + 1*CC_full, Y_resvar) +
        rep(alpha0[cluster_sizes > 20], cluster_sizes[cluster_sizes > 20])
      
      Y1 <- rep(NA, n)
      #Y1[sizecheck == 0] <- rnorm(sum(1-sizecheck), 0.5, Y_resvar) +
        #rep(alpha0[cluster_sizes <= 50], cluster_sizes[cluster_sizes <= 50])
      #Y1[sizecheck == 1] <- rnorm(sum(sizecheck), 1.5 + 2.25*(k-5), Y_resvar) +
        #rep(alpha0[cluster_sizes > 50], cluster_sizes[cluster_sizes > 50])
      sizeratio <- mean(cluster_sizes_full > 20) / mean(cluster_sizes_full <= 20)
      Y1[sizecheck == 0] <- rnorm(sum(1-sizecheck), 0 + 1*CC_full - 0.1*(k-1)*sizeratio,
                                  Y_resvar) +
        rep(alpha0[cluster_sizes <= 20], cluster_sizes[cluster_sizes <= 20])
      Y1[sizecheck == 1] <- rnorm(sum(sizecheck), 1 + 1*CC_full + 0.1*(k-1), Y_resvar) +
        rep(alpha0[cluster_sizes > 20], cluster_sizes[cluster_sizes > 20])
      
      # Exact 1:1 treatment allocation
      trtprob <- 0.5
      trt <- sample(1:num_clusters, num_clusters*trtprob)
      Z <- rep(0, n)
      Z[clusternum %in% trt] <- 1
      
      # Observed outcome and true values
      Y <- Z*Y1 + (1-Z)*Y0
      pate_true[i, floor(j/30), k] <- mean(Y1 - Y0)
      cate_true[i, floor(j/30), k] <- sum((Y1 - Y0) / cluster_sizes_full) / num_clusters
      
      # Estimation for Analysts 1 and 2
      cs_cent_full <- cluster_sizes_full - mean(cluster_sizes_full)
      cs_cent_full <- 1*(cluster_sizes_full > 20) - mean(1*(cluster_sizes_full > 20))
      pate_mod1 <- geeglm(Y ~ Z + cs_cent_full + CC_full,
                          corstr = "independence", id = clusternum)
      pate_mod2 <- geeglm(Y ~ Z + cs_cent_full + CC_full,
                          corstr = "exchangeable", id = clusternum)
      
      # Bias
      pate_bias1[i, floor(j/30), k] <-
        summary(pate_mod1)$coefficients[2, 1] - pate_true[i, floor(j/30), k]
      pate_bias2[i, floor(j/30), k] <-
        summary(pate_mod2)$coefficients[2, 1] - pate_true[i, floor(j/30), k]
      
      # SE
      pate_se1[i, floor(j/30), k] <- summary(pate_mod1)$coefficients[2, 2]
      pate_se2[i, floor(j/30), k] <- summary(pate_mod2)$coefficients[2, 2]
      
      # P-values
      pate_pval1[i, floor(j/30), k] <-
        summary(pate_mod1)$coefficients[2, 4]
      pate_pval2[i, floor(j/30), k] <-
        summary(pate_mod2)$coefficients[2, 4]
      
      # Estimation for Analyst 3
      # Model assisted test
      Yvec <- tapply(Y, clusternum, mean)
      Zvec <- tapply(Z, clusternum, mean)
      w <- cluster_sizes / n - 1 / num_clusters
      wYvec <- w*Yvec*num_clusters
      cs_cent <- cluster_sizes - mean(cluster_sizes)
      cs_cent <- 1*(cluster_sizes > 20) - mean(1*(cluster_sizes > 20))
      
      diff_pval <- coeftest(lm(wYvec ~ Zvec*CC + Zvec*cs_cent),
                            vcov = vcovHC(lm(wYvec ~ Zvec*CC + Zvec*cs_cent),
                                          type = "HC0"))[2, 4]
      
      # Second-stage which should have issues and not be done in practice
      # after discussion with co-authors
      if (diff_pval < 0.05) {
        pate_bias3[i, floor(j/30), k] <- pate_bias1[i, floor(j/30), k]
        pate_se3[i, floor(j/30), k] <- pate_se1[i, floor(j/30), k]
        pate_pval3[i, floor(j/30), k] <- pate_pval1[i, floor(j/30), k]
      } else if (diff_pval >= 0.05) {
        pate_bias3[i, floor(j/30), k] <- pate_bias2[i, floor(j/30), k]
        pate_se3[i, floor(j/30), k] <- pate_se2[i, floor(j/30), k]
        pate_pval3[i, floor(j/30), k] <- pate_pval2[i, floor(j/30), k]
      }
      
    }
    
  }
  
}

colMeans(pate_pval1[, 3, ] < 0.05)
colMeans(pate_pval3[, 3, ] < 0.05)
plot(1:9, colMeans(pate_pval1[, 3, ] < 0.05), ylim = c(0, 0.1))
lines(1:9, colMeans(pate_pval3[, 3, ] < 0.05))

# Make results figure
plotdat <- data.frame("Analysis" = c(rep("Independence GEE", 9),
                                    rep("Two-stage Procedure", 9)),
                      "Effect_Size" = rep(round(colMeans(cate_true[, 3, ]), 3), 2),
                      "T1E" = c(colMeans(pate_pval1[, 3, ] < 0.05),
                                 colMeans(pate_pval3[, 3, ] < 0.05)))

p1 <- ggplot(dat = plotdat[plotdat$Effect_Size > -0.35, ], aes(x = abs(Effect_Size),
                                                               y = T1E,
                          shape = Analysis, color = Analysis)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05) +
  geom_hline(yintercept = 0.05 + 1.96*sqrt(0.05*0.95/nsims), lty = 2) +
  geom_hline(yintercept = 0.05 - 1.96*sqrt(0.05*0.95/nsims), lty = 2) +
  scale_color_manual(values = c("Independence GEE" = "blue",
                                "Two-stage Procedure" = "red")) +
  xlab("Absolute c-ATE") +
  ylab("Type I Error for i-ATE test") + ylim(c(0, 0.1)) +
  theme_bw()

# Save results
save(pate_true, cate_true, pate_bias1, pate_bias2, pate_bias3,
     pate_se1, pate_se2, pate_se3, pate_pval1, pate_pval2, pate_pval3,
     file = "~/OneDrive_VUMC/Projects/ICS/Output/sims-scenario3.RData")


# Make results figure
plotdat2 <- data.frame("Analysis" = c(rep("Independence GEE", 9),
                                     rep("Two-stage Procedure", 9)),
                       "Effect_Size" = rep(round(colMeans(cate_true[, 3, ]), 3), 2),
                       "Bias" = c(colMeans(pate_bias1[, 3, ]),
                                  colMeans(pate_bias3[, 3, ])),
                       "MCSE" = c(colSds(pate_bias1[, 3, ]) / sqrt(nsims),
                                   colSds(pate_bias3[, 3, ]) / sqrt(nsims)),
                       "SE" = c(colMeans(pate_se1[, 3, ]),
                                colMeans(pate_se3[, 3, ])))
plotdat2$Lower <- plotdat2$Bias - 1.96* plotdat2$MCSE
plotdat2$Upper <- plotdat2$Bias + 1.96* plotdat2$MCSE

p2 <- ggplot(dat = plotdat2, aes(x = abs(Effect_Size), y = Bias,
                          shape = Analysis, color = Analysis)) +
  geom_line() +
  geom_point() +
  #facet_wrap(~Num_Clusters) +
  geom_hline(yintercept = 0) +
  #geom_hline(yintercept = 0.05 + 1.96*sqrt(0.05*0.95/nsims), lty = 2) +
  #geom_hline(yintercept = 0.05 - 1.96*sqrt(0.05*0.95/nsims), lty = 2) +
  # Add error shaded region with same colors as lines
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Analysis), alpha = 0.2) +
  scale_color_manual(values = c("Independence GEE" = "blue",
                                 "Two-stage Procedure" = "red")) +
  scale_fill_manual(values = c("Independence GEE" = "blue",
                               "Two-stage Procedure" = "red")) +
  xlab("Absolute c-ATE") +
  ylab("Bias for i-ATE") +
  theme_bw()

ggarrange(p1, p2, ncol = 1, nrow = 2)

# Relative efficiency table
se1_20 <- plotdat$SE[plotdat$Analyst == "Analyst 1" &
                     plotdat$Num_Clusters == 20][5]
se1_60 <- plotdat$SE[plotdat$Analyst == "Analyst 1" &
                     plotdat$Num_Clusters == 60][5]
se1_100 <- plotdat$SE[plotdat$Analyst == "Analyst 1" &
                      plotdat$Num_Clusters == 100][5]

se2_20 <- plotdat$SE[plotdat$Analyst == "Analyst 2" &
                     plotdat$Num_Clusters == 20][5]
se2_60 <- plotdat$SE[plotdat$Analyst == "Analyst 2" &
                     plotdat$Num_Clusters == 60][5]
se2_100 <- plotdat$SE[plotdat$Analyst == "Analyst 2" &
                      plotdat$Num_Clusters == 100][5]

se3_20 <- plotdat$SE[plotdat$Analyst == "Analyst 3" &
                     plotdat$Num_Clusters == 20][5]
se3_60 <- plotdat$SE[plotdat$Analyst == "Analyst 3" &
                     plotdat$Num_Clusters == 60][5]
se3_100 <- plotdat$SE[plotdat$Analyst == "Analyst 3" &
                      plotdat$Num_Clusters == 100][5]

re_2vs1_20 <- se2_20^2 / se1_20^2
re_2vs1_60 <- se2_60^2 / se1_60^2
re_2vs1_100 <- se2_100^2 / se1_100^2
re_3vs1_20 <- se3_20^2 / se1_20^2
re_3vs1_60 <- se3_60^2 / se1_60^2
re_3vs1_100 <- se3_100^2 / se1_100^2



plotdat2 <- cbind(pate_bias1[, 2, 4], pate_bias2[, 2, 4], pate_bias3[, 2, 4])
boxplot(plotdat2, beside=T)
