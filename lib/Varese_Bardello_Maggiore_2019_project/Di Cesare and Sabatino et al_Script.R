##### load data #####

ARGst <- read.csv("data_ARGsubtypes.csv", header=T, as.is=F)
ARGst<-ARGst[,-1]
MRGst <- read.csv("data_MRGsubtypes.csv", header=T, as.is=F)
MRGst<-MRGst[,-1]
ARGt <- read.csv("data_ARGtypes.csv", header=T, as.is=F)
MRGt <- read.csv("data_MRGtypes.csv", header=T, as.is=F)

##### load libraries #####

library(car)
library(emmeans)
library(ggplot2)
library(MASS)
library(performance)
library(vegan)
library(reshape2)
library("cowplot")

##### analyse data #####

# richness of ARGst

ARGst$samples <- factor(ARGst$samples, levels=c("LV", "RB1", "RB2", "RB3", "LM1", "LM2", "LM3"))
ARGst$month <- factor(ARGst$month, levels=c("August", "October", "December"))
ARGst$richness <- apply(ARGst[,3:length(names(ARGst))],1,function(x) sum(x>0))

ARGst_m <- MASS::glm.nb(richness ~ samples + month, data=ARGst)
performance::check_model(ARGst_m)
car::Anova(ARGst_m)
emmeans::emmeans(ARGst_m, pairwise ~ samples)->ph
ph
out_m<-capture.output(Anova(ARGst_m))
out_ph<-capture.output(ph$contrasts)
write(c("ARGs",out_m," ",out_ph), "stats.txt", append=T)

colS2<-c("#bcf60c","thistle3","darkslateblue")
g1<-ggplot2::ggplot(ARGst, aes(x=factor(samples,levels=c("LV", "RB1", "RB2", "RB3", "LM1", "LM2", "LM3")), y=richness))+
	geom_boxplot(lwd=0.5, color="black", alpha=0.5) + 
	geom_point(aes(fill=month, color=month),size=5, alpha=0.5)  +
	geom_line(aes(group=month,color=month), linetype=2) +
	xlab("") + ylab("Richness (antibiotic resistome)") + 
	scale_fill_manual(values=colS2)+
	scale_color_manual(values=colS2)+
	theme(legend.position = c(0.85, 0.80))
g1

# differences in ARGst composition

betabrayA<-vegdist(ARGst[,c(3:810)],method="bray")
plot(hclust(betabrayA, method="average"),hang=-1,  sub='', xlab='', cex=0.5) #plot cluster analysis of betapair

adonis<-adonis(betabrayA~ARGst$samples+ARGst$month, permutations=9999)
adonis
out_a<-capture.output(adonis)
write(out_a, "stats.txt", append=T)

NMDS<-metaMDS(ARGst[,c(3:810)], distance="bray", k=3)
plot(NMDS)

colNew2<-c("thistle3",'#bcf60c','#3cb44b','#808000','steelblue1',"cadetblue4", "darkslateblue")
data.scores <- as.data.frame(scores(NMDS),display="sites")  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores

veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame() #sets up a data frame before running the function.
for(g in levels(data.scores$Site)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(data.scores [data.scores$Site==g,],
                                                         veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,Site=g))
}

Site<-factor(ARGst$samples,levels=c("LV", "RB1", "RB2", "RB3", "LM1","LM2", "LM3"))
cmp1<-ggplot() + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,color=Site,fill=Site),size=5) + # add the point markers
  scale_fill_manual(values=colNew2) +
  scale_color_manual(values=colNew2)+
  scale_shape_manual(values=c(23,21,24,22))+ 
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,color=Site),alpha=0.4, size=0.5, linetype=1)+
  theme_bw()
cmp1

# abundance of ARGst

ARGst$samples <- factor(ARGst$samples,levels=c("LV", "RB1", "RB2", "RB3", "LM1","LM2", "LM3"))
ARGst$month <- factor(ARGst$month, levels=c("August", "October", "December"))
ARGst$abundance <- rowSums(ARGst[,3:(length(names(ARGst))-1)])

ARGst_a <- lm(abundance ~ samples + month, data=ARGst)
performance::check_model(ARGst_a)
car::Anova(ARGst_a)
emmeans::emmeans(ARGst_a, pairwise ~ samples)->phA
emmeans::emmeans(ARGst_a, pairwise ~ month)->phB
phA
phB
out_a<-capture.output(Anova(ARGst_a))
out_phA<-capture.output(phA$contrasts)
out_phB<-capture.output(phB$contrasts)
write(c(out_a," ",out_phA," ",out_phB), "stats.txt", append=T)

ggplot2::ggplot(ARGst, aes(samples, abundance)) +
	geom_boxplot(lwd=0.5, color="light grey") +
	xlab("") + ylab("") + ggtitle("ARGst") + 
	geom_point(aes(fill=month, group=month), size=4, shape=21) +
	geom_line(aes(group=month,color=month), linetype=2) +
	scale_fill_manual(values=c("#F99E1F", "#990000", "#040273")) +
	scale_color_manual(values=c("#F99E1F", "#990000", "#040273")) + 
	theme_classic()


# richness of MRGst

MRGst$samples <- factor(MRGst$samples, levels=c("LV", "RB1", "RB2", "RB3", "LM1","LM2", "LM3"))
MRGst$month <- factor(MRGst$month, levels=c("August", "October", "December"))
MRGst$richness <- apply(MRGst[,3:length(names(MRGst))],1,function(x) sum(x>0))

MRGst_m <- MASS::glm.nb(richness ~ samples + month, data=MRGst)
performance::check_model(MRGst_m)
car::Anova(MRGst_m)
emmeans::emmeans(MRGst_m, pairwise ~ samples)->ph
ph
out_m<-capture.output(Anova(MRGst_m))
out_ph<-capture.output(ph$contrasts)
write(c(" ","MRGs",out_m," ",out_ph), "stats.txt", append=T)

g2<-ggplot2::ggplot(MRGst, aes(x=factor(samples,levels=c("LV", "RB1", "RB2", "RB3", "LM1", "LM2", "LM3")), y=richness))+
	geom_boxplot(lwd=0.5, color="black", alpha=0.5) + 
	geom_point(aes(fill=month, color=month),size=5, alpha=0.5)  +
	geom_line(aes(group=month,color=month), linetype=2) +
	xlab("") + ylab("Richness (metal resistome)") + 
	scale_fill_manual(values=colS2)+
	scale_color_manual(values=colS2)+
	theme(legend.position = "none")
g2

plot_grid(g1, g2, labels=c("A","B"), align="v",ncol=1)

# differences in MRGst composition

betabrayM<-vegdist(MRGst[,c(3:278)],method="bray")
plot(hclust(betabrayM, method="average"),hang=-1,  sub='', xlab='', cex=0.5) #plot cluster analysis of betapair

adonis<-adonis(betabrayM~MRGst$samples+MRGst$month, permutations=9999)
adonis
out_a<-capture.output(adonis)
write(out_a, "stats.txt", append=T)

NMDS2<-metaMDS(MRGst[,c(3:278)], distance="bray", k=3)
plot(NMDS2)
data.scores2 <- as.data.frame(scores(NMDS2),display="sites")  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores2$site <- rownames(data.scores2)  # create a column of site names, from the rownames of data.scores

df_ell2 <- data.frame() #sets up a data frame before running the function.
for(g in levels(data.scores2$Site)){
  df_ell2 <- rbind(df_ell2, cbind(as.data.frame(with(data.scores2 [data.scores2$Site==g,],
                                                         veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,Site=g))
}

Site<-factor(MRGst$samples,levels=c("LV", "RB1", "RB2", "RB3", "LM1","LM2", "LM3"))
cmp2<-ggplot() + 
  geom_point(data=data.scores2,aes(x=NMDS1,y=NMDS2,color=Site,fill=Site),size=5) + # add the point markers
  scale_fill_manual(values=colNew2) +
  scale_color_manual(values=colNew2)+
  scale_shape_manual(values=c(23,21,24,22))+
  geom_path(data=df_ell2, aes(x=NMDS1, y=NMDS2,color=Site),alpha=0.4, size=0.5, linetype=1)+
  theme_bw()
cmp2

plot_grid(cmp1, cmp2, labels=c("A","B"), align="v",ncol=1)

# abundance of MRGst

MRGst$samples <- factor(MRGst$samples, levels=c("LV", "RB1", "RB2", "RB3", "LM1","LM2", "LM3"))
MRGst$month <- factor(MRGst$month, levels=c("August", "October", "December"))
MRGst$abundance <- rowSums(MRGst[,3:(length(names(MRGst))-1)])

MRGst_a <- lm(abundance ~ samples + month, data=MRGst)
performance::check_model(MRGst_a)
car::Anova(MRGst_a)
emmeans::emmeans(MRGst_a, pairwise ~ samples)->phA
emmeans::emmeans(MRGst_a, pairwise ~ month)->phB
phA
phB
out_a<-capture.output(Anova(MRGst_a))
out_phA<-capture.output(phA$contrasts)
out_phB<-capture.output(phB$contrasts)
write(c(out_a," ",out_phA," ",out_phB), "stats.txt", append=T)

ggplot2::ggplot(MRGst, aes(samples, abundance)) +
	geom_boxplot(lwd=0.5, color="light grey") +
	xlab("") + ylab("") + ggtitle("MRGst") + 
	geom_point(aes(fill=month, group=month), size=4, shape=21) +
	geom_line(aes(group=month,color=month), linetype=2) +
	scale_fill_manual(values=c("#F99E1F", "#990000", "#040273")) +
	scale_color_manual(values=c("#F99E1F", "#990000", "#040273")) + 
	theme_classic()

#  correlation

MRGst <- read.csv("data_MRGsubtypesB.csv", header=T, as.is=F)
MRGst<-MRGst[,-1]
MRGst$samples <- factor(MRGst$samples, levels=c("LV", "RB1", "RB2", "RB3", "LM1","LM2", "LM3"))
MRGst$month <- factor(MRGst$month, levels=c("August", "October", "December"))
MRGst$abundance <- rowSums(MRGst[,3:(length(names(MRGst))-1)])
pearson<-cor.test(ARGst$abundance,MRGst$abundance,method ="pearson")
pearson
plot(ARGst$abundance~MRGst$abundance)

betabrayM<-vegdist(MRGst[,c(3:264)],method="bray")
mtl<-mantel(betabrayA,betabrayM,method="pearson")
mtl

##########  ARG and MRG type
ARGt$samples <- factor(ARGt$samples, levels=c("LV", "RB1", "RB2", "RB3", "LM1","LM2", "LM3"))
ARGt$month <- factor(ARGt$month, levels=c("August", "October", "December"))
rownames(ARGt)<-paste(ARGt$samples,"",ARGt$month)
MRGt$samples <- factor(MRGt$samples, levels=c("LV", "RB1", "RB2", "RB3", "LM1","LM2", "LM3"))
MRGt$month <- factor(MRGt$month, levels=c("August", "October", "December"))
rownames(MRGt)<-paste(MRGt$samples,"",MRGt$month)

AC<-ARGt[,c(3:28)]
tAC<-t(AC)
View(tAC)
MC<-MRGt[,c(3:20)]
tMC<-t(MC)
View(tMC)
allST1<-as.data.frame(tAC)
datm1<-melt(cbind(allST1, ind=row.names(allST1), id.vars=c('ind')))

colNew3<-c("darkslateblue",	'#4363d8',"slategray2","slategray4", 'grey83',"snow4",'steelblue1','#e6194b',
	'#800000','#f58231','#fffac8',"red",'#ffe119',"pink3",'#fabebe','#ffd8b1',"peachpuff3",	'#aaffc3',			
	'#bcf60c',	'#3cb44b',	'#808000',"brown","cadetblue4","thistle3","mediumorchid",'#e6beff')
a<-ggplot(datm1,aes(x = factor(variable, level=c("LV  August","RB1  August","RB2  August","RB3  August","LM1  August","LM2  August","LM3  August","LV  October","RB1  October","RB2  October",
	"RB3  October","LM1  October","LM2  October","LM3  October","LV  December","RB1  December","RB2  December","RB3  December","LM1  December","LM2  December","LM3  December")), 
	y = value,fill = ind)) + geom_bar(position="stack", stat="identity") +scale_fill_manual(values= colNew3)+ guides(fill=guide_legend(ncol=7))+theme(text = element_text(size=15), 
	axis.text.x = element_text(angle=90, vjust=0), legend.position="bottom", legend.text=element_text(size=10), legend.title=element_blank()) + scale_x_discrete(name="") + 
	scale_y_continuous(name="Relative Abundance") + theme(axis.title.y=element_text(colour="#999999"))+labs(tag ="A")
a

allST2<-as.data.frame(tMC)
datm2<-melt(cbind(allST2, ind=row.names(allST2), id.vars=c('ind')))
b<-ggplot(datm2,aes(x = factor(variable, level=c("LV  August","RB1  August","RB2  August","RB3  August","LM1  August","LM2  August","LM3  August","LV  October","RB1  October","RB2  October",
	"RB3  October","LM1  October","LM2  October","LM3  October","LV  December","RB1  December","RB2  December","RB3  December","LM1  December","LM2  December","LM3  December")),  
	y = value,fill = ind)) + geom_bar(position="stack",stat="identity") +scale_fill_manual(values= colNew)+ guides(fill=guide_legend(ncol=6))+theme(text = element_text(size=15), 
	axis.text.x = element_text(angle=90, vjust=0), legend.position="bottom", legend.text=element_text(size=10), legend.title=element_blank()) + scale_x_discrete(name="") + 
	scale_y_continuous(name="Relative Abundance") + theme(axis.title.y=element_text(colour="#999999"))+labs(tag ="B")
b

plot_grid(a, b, align="v",ncol=1)