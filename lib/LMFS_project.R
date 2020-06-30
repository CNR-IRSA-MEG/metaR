

write("Info for article:", "info.csv") 
version<-R.Version() 
write(version$version.string, "info.csv", append=T) 
getwd()
###Prepare OTU table and taxonomy file

OTUListM<-read.delim("datosM.csv", sep=";")
dim(OTUListM) 
names(OTUListM)
rownames(OTUListM)<-OTUListM[,1] 
OTUListM<-OTUListM[,-1] 

###Normalise read numbers:
readnr<- colSums(OTUListM) 
colSums(OTUListM)
summary(readnr) 
hist(readnr)
samples<-colnames(OTUListM)
sampleNr<-cbind(samples,readnr) 
sampleNr
sampleNr<-cbind(colnames(OTUListM),colSums(OTUListM)) #the same thing as the three rows above but coded in a more effient way.
 


tOTUListM<-as.data.frame(t(OTUListM)) #Transpose 
dim(tOTUListM)

#Read in the variables
variables<-read.delim("variables.csv",sep=";")


#######################################

write("statistical results:", "statsM.csv")

### ALPHA Diversity

library(vegan)
richness<-specnumber(tOTUListM) 
shannon<-diversity(tOTUListM)
evenness<- shannon/log(specnumber(tOTUListM))


Alpha<-cbind(richness,shannon,evenness)
write.csv(Alpha, "alphaM.csv") 
str(variables)
glmRichM<-glm(richness~variables$matrix, family=poisson(link="log"))
summary(glmRichM)
glmRichOUTM<-capture.output(summary(glmRichM)) 

write("differences in richness M matrix:", "stats.csv", append=T)
write(glmRichOUTM, "statsM.csv", append=T)


###BETADIVERSITY:

beta<-vegdist(tOTUListM, method="bray", diag=T, upper=T) 
beta

plot(hclust(beta, method="average"), main="bray curtis", hang=-1, cex=0.7)

adonis<-adonis(beta~variables$matrix, strata = variables$replicate)
adonis
adonisOUT<-capture.output(adonis)
write.csv(adonisOUT, "adonisM.csv")

library(RVAideMemoire)
Matrix<-pairwise.perm.manova(beta, variables$matrix)
Matrix
write("pairwise matrix Mbetadiversity", "stats.csv", append=T)
write(Matrix$p.value, "stats.csv", append=T)


library("reshape2")

#Plot as heat map
HEATM<-read.delim("heatM.csv", sep=";") 
dim(HEATM) 
names(HEATM)
rownames(HEATM)<-HEATM[,1]
HEATM<-HEATM[,-1]
mST<-as.matrix(HEATM)
library(RColorBrewer)

heatmap(mST,col=colorRampPalette(brewer.pal(9, "Blues"))(9), scale="row", Rowv = NA)

rangeLegendUnit <- round((max(0.659)-min(0))/9, digits=3)
Legend1 <- paste("0-", rangeLegendUnit)
Legend2 <- paste(rangeLegendUnit, "-", rangeLegendUnit*2)
Legend3 <- paste(rangeLegendUnit*2, "-", rangeLegendUnit*3)
Legend4 <- paste(rangeLegendUnit*3, "-", rangeLegendUnit*4)
Legend5 <- paste(rangeLegendUnit*4, "-", rangeLegendUnit*5)
Legend6 <- paste(rangeLegendUnit*5, "-", rangeLegendUnit*6)
Legend7 <- paste(rangeLegendUnit*6, "-", rangeLegendUnit*7)
Legend8 <- paste(rangeLegendUnit*7, "-", rangeLegendUnit*8)
Legend9 <- paste(rangeLegendUnit*8, "-", rangeLegendUnit*9)

legend(x="topleft", cex =0.7, legend=c(Legend1, Legend2, Legend3, Legend4, Legend5, Legend6, Legend7, Legend8, Legend9),
       fill=colorRampPalette(brewer.pal(9, "Blues"))(9))

######ANOVA

library(stats)     
library(agricolae)  
library(gplots)     

BASE <- read.csv(file = "anovaM.csv", sep = ";", header = TRUE)
BASE

NUM_COL_VR <- 1
NUM_COL_VR

NUM_COL_FACTOR <- 2
NUM_COL_FACTOR
NUM_COL_SEL <- c(NUM_COL_VR , NUM_COL_FACTOR)
NUM_COL_SEL
VR <- BASE[,1]
FACTOR <- BASE[,2]
VR
FACTOR

is.factor(FACTOR)
is.numeric(VR)

FACTOR <- as.factor(FACTOR)

LM <- lm(VR ~ FACTOR)
LM

ANOVA_COMPLETO <- aov(LM)
ANOVA_COMPLETO

TABLA_ANOVA <- summary(ANOVA_COMPLETO)[[1]]
TABLA_ANOVA

TUKEY_COMPLETO <- HSD.test(LM,"FACTOR")

TABLA_TUKEY <- TUKEY_COMPLETO$groups
TABLA_TUKEY

plot(TUKEY_COMPLETO)

plotmeans(VR ~ FACTOR, col="blue")

RESIDUOS <- ANOVA_COMPLETO$residuals
RESIDUOS
SHAPIRO <- shapiro.test(RESIDUOS)
SHAPIRO
BARTLETT <- bartlett.test(RESIDUOS, FACTOR)
BARTLETT

############   FIN DEL SCRIPT  ############

