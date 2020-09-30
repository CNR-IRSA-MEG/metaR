added <- read.csv("minor dataset.csv", as.is=F)
names(added)
attach(added)

#################
##### plots #####
#################

par(mai=c(0,1,0.5,0), cex.axis=1.2, cex.lab=1.5)

i <- sequences
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(15,max(i)), pch=21, bg=as.factor(origin), cex=sequences/12, ylab="number of sequences")

i <- haplotypes
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(5,30), pch=21, bg=as.factor(origin), cex=sequences/12, ylab="number of haplotypes")

i <- haplotypes/sequences
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0.2,1), pch=21, bg=as.factor(origin), cex=sequences/12, ylab="number of haplotypes / number of sequences")

i <- median.genetic.distance
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0,0.2), pch=21, bg=as.factor(origin), cex=sequences/12, ylab="median genetic distances")

i <- GMYC
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(0,12), pch=21, bg=as.factor(origin), cex=sequences/12, ylab="number of GMYC units")

i <- Tajima.D
stripchart(i, method="jitter", jitter = 0.5, vertical=TRUE, frame=F, ylim=c(-2,1), pch=21, bg=as.factor(origin), cex=sequences/12, ylab="Tajima's D")

##################
##### models #####
##################

library(rsq)

i <- sequences
model_GLMp <- glm(log(i) ~ origin, family=quasipoisson)
summary(model_GLMp)
rsq.partial(model_GLMp, adj=T)

i <- haplotypes
model_GLMp <- glm(log(i) ~ origin + sequences, family=quasipoisson)
summary(model_GLMp)
rsq.partial(model_GLMp, adj=T)

model_GLMb <- glm(cbind(haplotypes, sequences) ~ origin + sequences, family=binomial)
summary(model_GLMb)
rsq.partial(model_GLMb, adj=T)

i <- median.genetic.distance
model_LM <- lm(i ~ origin + sequences)
summary(model_LM)
rsq.partial(model_LM, adj=T)

i <- GMYC
model_GLMp <- glm(log(i+1) ~ origin + sequences, family=quasipoisson)
summary(model_GLMp)
rsq.partial(model_GLMp, adj=T)

i <- Tajima.D
model_LM <- lm(i ~ origin + sequences)
summary(model_LM)
rsq.partial(model_LM, adj=T)