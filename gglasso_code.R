# group lasso regression with stability selection #
# example with M1-CB Filigree subnetworks from MMIP data #
# same code was used for M3-CB and M1-M3 #


setwd("")

source('lib2/stabsel_gglasso.R')

library(gglasso)
library(ggplot2)

meta.dat <- read.csv("MMIP covariates.csv", 
                     header = T, check.names = F, stringsAsFactors = F)
meta.dat <- meta.dat[,c(1,2,6,11,18)]
meta.dat <- meta.dat[-c(which(is.na(meta.dat$BMI_baseline)), which(is.na(meta.dat$Fenton_BW_Z))),]

m1.cb <- read.csv("Age_BMI_adjusted_M1_CB_final_294_collpased_scaled.csv", 
                  header = T, check.names = F, stringsAsFactors = F, row.names = 1)

nodelist <- read.csv("Age_BMI_adjusted_M1_CB_final_294_collpased_scaled_nodelist_nrep50.csv", 
                     header = T, stringsAsFactors = F, check.names = F)
nodelist <- nodelist[,c(1,10)]

m1.dat <- m1.cb[m1.cb$Group == "M1",-1]
cb.dat <- m1.cb[m1.cb$Group == "CB",-1]

m1.clusters <- merge(nodelist, t(m1.dat), by.x = "Compound", by.y = "row.names", sort = F)
rownames(m1.clusters) <- m1.clusters[,1]
m1.clusters <- m1.clusters[,-1]
rownames(m1.clusters) <- paste0(rownames(m1.clusters), " (M1)")
m1.clusters$membership[which(is.na(m1.clusters$membership))] <- 21
colnames(m1.clusters)[-1] <- substr(colnames(m1.clusters)[-1], 4, 
                                    nchar(colnames(m1.clusters)[-1]))

cb.clusters <- merge(nodelist, t(cb.dat), by.x = "Compound", by.y = "row.names", sort = F)
rownames(cb.clusters) <- cb.clusters[,1]
cb.clusters <- cb.clusters[,-1]
rownames(cb.clusters) <- paste0(rownames(cb.clusters), " (CB)")
cb.clusters$membership[which(is.na(cb.clusters$membership))] <- 21
cb.clusters$membership <- cb.clusters$membership + max(m1.clusters$membership)
colnames(cb.clusters)[-1] <- substr(colnames(cb.clusters)[-1], 4, 
                                    nchar(colnames(cb.clusters)[-1]))

m1.cb.clusters <- rbind.data.frame(m1.clusters, cb.clusters)
m1.cb.clusters <- m1.cb.clusters[,(colnames(m1.cb.clusters)[-1] %in% meta.dat$ID)]
m1.cb.cluster.num <- m1.cb.clusters$membership
m1.cb.clusters.for.gglasso <- t(m1.cb.clusters[,-1])

################################################
###########################################

stabsel_m1_cb_gglasso <- stabsel_gglasso(X = as.matrix(m1.cb.clusters.for.gglasso),
					Y = meta.dat$Fenton_BW_Z,
					nSample = 500,
					group_membership = m1.cb.cluster.num,
					method = "ls",
					seed.add = 1)

stabsel_m1_cb_lambda_rank <- matrix(0, 42, 1000)
for (i in 1:ncol(stabsel_m1_cb_lambda_rank)){
	stabsel_m1_cb_lambda_rank[,i] <- as.numeric(factor(rank(-stabsel_m1_cb_gglasso[,i])))
}
rank_freq <- matrix(0, nrow = 42, ncol = 5)
for ( i in 1:nrow(stabsel_m1_cb_lambda_rank)){
	rank_freq[i,1] <- length(which(stabsel_m1_cb_lambda_rank[i,] == 1))
	rank_freq[i,2] <- length(which(stabsel_m1_cb_lambda_rank[i,] == 2))
	rank_freq[i,3] <- length(which(stabsel_m1_cb_lambda_rank[i,] == 3))
	rank_freq[i,4] <- length(which(stabsel_m1_cb_lambda_rank[i,] == 4))
	rank_freq[i,5] <- length(which(stabsel_m1_cb_lambda_rank[i,] == 5))
}
rownames(rank_freq)[1:21] <- paste0("Cluster_", 1:21, " (M1)")
rownames(rank_freq)[22:42] <- paste0("Cluster_", 1:21, " (CB)")
rownames(rank_freq)[1:21] <- paste0("Cluster_", 1:21, " (M1)")

rank_freq <- cbind.data.frame(rank_freq, as.numeric(lapply(1:42, function (v) mean(rank_freq[v,]))))
colnames(rank_freq)[6] <- "mean_rank"
rank_freq <- rank_freq[rev(order(rank_freq$mean_rank)),]
rank_freq2 <- cbind.data.frame(rownames(rank_freq), rank_freq$Rank_1)
colnames(rank_freq2) <- c("Subnetwork #", "Rank_1_frequency")
colnames(rank_freq)[1:5] <- paste0("Rank_", 1:5)
rank_freq2 <- cbind.data.frame(rownames(rank_freq), rank_freq$Rank_1)
colnames(rank_freq2) <- c("Subnetwork #", "Rank_1_frequency")

#################################################################################

png("gglasso_stabsel_M1_CB_glassoOnly_nReps50.png", 
    height = 7, width = 10, units = "in", res = 300)

ggplot(data = rank_freq2, aes(x = reorder("`Subnetwork #`", -Rank_1_frequency), 
                              y = Rank_1_frequency)) + 
  geom_bar(stat = "identity", color = "black", fill ="#FF6666") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title.x = element_blank(), legend.position = "none")

dev.off()


write.csv(rank_freq2, file = "rank_freq2.csv")
