anta <- read.csv("anta data.csv", as.is=F, header=T)
names(anta)
summary(anta)

attach(anta)

#################
##### plots #####
#################

par(mai=c(0,1,0.5,0), cex.axis=1.5, cex.lab=2)

i <- geoD_max
stripchart(i, method="jitter", jitter = 0.1, vertical=TRUE, frame=F, ylim=c(0,max(i, na.rm=T)), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="maximum geographic distance (km)")

i <- N_sequences
stripchart(i, method="jitter", jitter = 10, vertical=TRUE, frame=F, ylim=c(0,max(i)), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="number of sequences")

i <- N_haplotypes
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0,max(i)), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="number of haplotypes")

i <- N_haplotypes/N_sequences
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0,1), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="haplotype / sequence number")

i <- N_ABGD
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0,max(i)), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="number of ABGD units")

i <- min_derivative
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0,max(i)), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="minimum derivative of accumulation curve")

i <- N_GMYC
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0,max(i)), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="number of GMYC units")

i <- genD_median
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0,0.2), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="median genetic distances")

i <- genD_mean
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0,max(i)), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="mean genetic distances")

i <- genD_max
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0,max(i)), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="maximum genetic distances")

i <- MRCA
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(-1,-0.4), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="relative crown age")

i <- geoD_median
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0,2000), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="median geographic distances")

i <- geoD_mean
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0,max(i)), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="mean geographic distances")

i <- geoD_max
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0,max(i)), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="maximum geographic distances")

i <- Mantel_r
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0,0.7), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="r value")

i <- Tajima_D
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(-4,0), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="Tajima's D")

i <- Fu_Li_F
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(-1.5,1.5), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="Fu and Li's F")

i <- Fst
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0,0.6), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="Fst")

i <- ANeCA_negatives/(ANeCA_negatives+ANeCA_positives)
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0.6,1), pch=21, bg=as.factor(origin), cex=log(N_sequences)-2, ylab="proportion of ANeCA non inferred steps")

##################
##### models #####
##################

names(anta)
library(psych)
anta_corr <- anta[,c(4,7,25,9:17,19,21,23)]
colnames(anta_corr) <- c("Sequences", "Haplotypes", "First derivative", "ABGD", "GMYC", "Median genetic", "Max genetic", "MRCA", "Median geographic", "Mean geographic", "Max geographic", "Mantel r", "Tajima's D", "Fu & Li's F", "Fst")
pairs.panels(anta_corr)

library(rsq)

i <- N_sequences
model_GLMp <- glm(log(i) ~ origin, family=quasipoisson)
summary(model_GLMp)
rsq.partial(model_GLMp, adj=T)

i <- N_haplotypes
model_GLMp <- glm(log(i) ~ origin + N_sequences, family=quasipoisson)
summary(model_GLMp)
rsq.partial(model_GLMp, adj=T)

model_GLMb <- glm(cbind(N_haplotypes, N_sequences) ~ origin + N_sequences, family=binomial)
summary(model_GLMb)
rsq.partial(model_GLMb, adj=T)

i <- min_derivative
model_LM <- lm(i ~ origin + N_sequences)
summary(model_LM)
rsq.partial(model_LM, adj=T)

i <- N_ABGD
model_GLMp <- glm(log(i+1) ~ origin + N_sequences, family=quasipoisson)
summary(model_GLMp)
rsq.partial(model_GLMp, adj=T)

i <- N_GMYC
model_GLMp <- glm(log(i+1) ~ origin + N_sequences, family=quasipoisson)
summary(model_GLMp)
rsq.partial(model_GLMp, adj=T)

i <- genD_median
model_LM <- lm(i ~ origin + N_sequences)
summary(model_LM)
rsq.partial(model_LM, adj=T)

i <- genD_max
model_LM <- lm(i ~ origin + N_sequences)
summary(model_LM)
rsq.partial(model_LM, adj=T)

i <- MRCA
model_LM <- lm(i ~ origin + N_sequences)
summary(model_LM)
rsq.partial(model_LM, adj=T)

i <- Tajima_D
model_LM <- lm(i ~ origin + N_sequences)
summary(model_LM)
rsq.partial(model_LM, adj=T)


##### subset = main

anta_main <- subset(anta, dataset=="main", drop=T)
summary(anta_main)
attach(anta_main)

i <- N_sequences
model_GLMp <- glm(log(i) ~ origin, family=quasipoisson)
summary(model_GLMp)
rsq.partial(model_GLMp, adj=T)

i <- geoD_max
model_LM <- lm(i ~ origin + N_sequences)
summary(model_LM)
rsq.partial(model_LM, adj=T)

i <- N_haplotypes
model_GLMp <- glm(log(i) ~ origin + N_sequences + geoD_max, family=quasipoisson)
summary(model_GLMp)
rsq.partial(model_GLMp, adj=T)

model_GLMb <- glm(cbind(N_haplotypes, N_sequences) ~ origin + N_sequences + geoD_max, family=binomial)
summary(model_GLMb)
rsq.partial(model_GLMb, adj=T)

i <- min_derivative
model_LM <- lm(i ~ origin + N_sequences + geoD_max)
summary(model_LM)
rsq.partial(model_LM, adj=T)

i <- N_ABGD
model_GLMp <- glm(log(i+1) ~ origin + N_sequences + geoD_max, family=quasipoisson)
summary(model_GLMp)
rsq.partial(model_GLMp, adj=T)

i <- N_GMYC
model_GLMp <- glm(log(i+1) ~ origin + N_sequences + geoD_max, family=quasipoisson)
summary(model_GLMp)
rsq.partial(model_GLMp, adj=T)

i <- genD_median
model_LM <- lm(i ~ origin + N_sequences + geoD_max)
summary(model_LM)
rsq.partial(model_LM, adj=T)

i <- genD_max
model_LM <- lm(i ~ origin + N_sequences + geoD_max)
summary(model_LM)
rsq.partial(model_LM, adj=T)

i <- MRCA
model_LM <- lm(i ~ origin + N_sequences + geoD_max)
summary(model_LM)
rsq.partial(model_LM, adj=T)

i <- Mantel_r
model_LM <- lm(i ~ origin + N_sequences + geoD_max)
summary(model_LM)
rsq.partial(model_LM, adj=T)

i <- Tajima_D
model_LM <- lm(i ~ origin + N_sequences + geoD_max)
summary(model_LM)
rsq.partial(model_LM, adj=T)

i <- Fu_Li_F
model_LM <- lm(i ~ origin + N_sequences + geoD_max)
summary(model_LM)
rsq.partial(model_LM, adj=T)

i <- Fst
model_LM <- lm(i ~ origin + N_sequences + geoD_max)
summary(model_LM)
rsq.partial(model_LM, adj=T)

model_GLMb <- glm(cbind(ANeCA_negatives, ANeCA_positives) ~ origin + N_sequences + geoD_max, family=binomial)
summary(model_GLMb)
rsq.partial(model_GLMb, adj=T)

