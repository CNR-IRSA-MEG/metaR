
setwd("C:/Users/Ester/Desktop/CNR/WWTP/Microplastic/") #define the directory where you keep you script and your data

write("Info for article:", "info.csv") #write some information to a file that you will use for the article.
version<-R.Version() #Fist information which R version are you using?
write(version$version.string, "info.csv", append=T) #Write it to the file.

###Prepare OTU table and taxonomy file

OTUList<-read.csv("otutable.csv") #reas in OTU table in csv formate for txt use read.delim
rownames(OTUList)<-OTUList[,1] #OTU names will become row names, [,1] means all the rows from the first col
OTUList<-OTUList[,-1] #remove OTU names form the dataframe

###Normalise read numbers:
readnr<- colSums(OTUList) #make the sum of each column to get the number of reads per sample
summary(readnr) 
samples<-colnames(OTUList) #extract the colnames from the dataframe to get a vector with the names of the samples
sampleNr<-cbind(samples,readnr) #make a new dataframe by binding the two vectors as columns 

sampleNr<-cbind(colnames(OTUList),colSums(OTUList)) #the same thing as the three rows above but coded in a more effient way.
  #this is infact your first result that you will need for publication
write.csv(sampleNr, "tableS1.csv") #Write the read numers to a file which will be saved in your working direcotry

tOTUList<-as.data.frame(t(OTUList)) #Transpose the OTU table! 

library("GUniFrac")
rareOTU<-Rarefy(tOTUList, depth = min(readnr))  #Rarefaction to number of reads of the sample with the least reads
rffOTU<-as.data.frame(rareOTU$otu.tab.rff) #the output is an array 

summary(rowSums(rffOTU)) #if you want to check if it worked make the summary per row they should all be the same...

min<-min(readnr)
write("rarefied to nr of reads:", "info.csv", append=TRUE)
write(min, "info.csv", append=TRUE)


raOTU <- rffOTU[ ,colSums(rffOTU)!=0] #some otus don't have any reads anymore, lets remove them: != means are not equal
traOTU<-t(raOTU) #transpose dataframe 
traOTU<-traOTU[order(row.names(traOTU)),] #sort by OTU number

###READ IN TAXONOMY FILE AND ADJUST OTU TABLE

taxonomy<-read.csv("taxonomy.csv") #taxonomy file 
row.names(taxonomy)<-taxonomy[,1] #set row names
taxonomy<-taxonomy[,-1]
STaxa<-taxonomy[order(row.names(taxonomy)),] #!!!Make sure that the OTUs in the table and in the taxonomy file are in the same order
taxaall<-STaxa[row.names(STaxa)%in%row.names(traOTU),] #We removed some OTUs now lets subset taxa to the ones that are still in the OTU table
taxaNC<-subset(taxaall, taxaall$c!="c:Chloroplast") #Remove OTUs from chloroplasts 16S
taxaNM<-subset(taxaall, taxaall$f!="f:Mitochondria") #Remove OTUs from mitochondial 16S
taxa<-subset(taxaNM, taxaNM$p!="") #Let's kick out OTUs that are not at least identified to a Phylum level
traOTU<-as.data.frame(traOTU[row.names(traOTU)%in%row.names(taxa),] ) #Keep OTUs that are in taxa
raOTU<-as.data.frame(t(traOTU)) #remake also the transposed dataframe
samples_otus<-dim(raOTU)
write("Nr of samples and Nr of otus left after cleaning:", "info.csv", append=TRUE)
write(samples_otus, "info.csv", append=TRUE)

#MAKE A PRESENCE ABSENCE TABLE
paOTU<-raOTU #make a new dataframe
paOTU[paOTU>0]=1 #in the new dataframe replace everything that it above 0 with 1

#Read in the variables
variables<-read.csv("variables.csv") #Read in your file with the variables per samples

#######################################
######START ANALYSIS OF DATA!
write("statistical results:", "stats.csv")

### ALPHA Diversity

library(vegan) #vegan has a very good tutorial pdf /documentation
richness<-specnumber(raOTU) #Number of species 
shannon<-diversity(raOTU)
evenness<- shannon/log(specnumber(raOTU))

par(mfrow=c(3,1)) #split the plot window in 3
plot(richness~variables$matrix, xlab="matrix")
plot(shannon~variables$matrix, xlab="matrix")
plot(evenness~variables$matrix, xlab="matrix")
dev.off()
Alpha<-cbind(richness,shannon,evenness)
write.csv(Alpha, "alpha.csv") #this is your second result

glmRich<-glm(richness~variables$matrix, family=poisson(link="log")) #Is there a statistical difference in richness between Microplastic and water samples?
summary(glmRich)
glmRichOUT<-capture.output(summary(glmRich)) #Write output in a way that is easily writable in a file

write("differences in richness MP vs water:", "stats.csv", append=T)
write(glmRichOUT, "stats.csv", append=T)


###BETADIVERSITY:

beta<-vegdist(raOTU, method="bray") #calculate beta diversit with bray curtis dissimilarity

plot(hclust(beta, method="average"), main="bray curtis", hang=-1, cex=0.7)

#PERMANOVA analysis: How much of the variance in betadiversity is explained by different variables
adonis<-adonis(beta~variables$matrix+variables$desin+variables$date, strata = variables$replicate) #if you want to considere interactions: variables$matrix*variables$desin*variables$date
adonis
adonisOUT<-capture.output(adonis)
write.csv(adonisOUT, "adonisAll.csv")

library(RVAideMemoire)
MPfilter<-pairwise.perm.manova(beta, variables$matrix)
MPfilter
write("pairwise MP filter betadiversity", "stats.csv", append=T)
write(MPfilter$p.value, "stats.csv", append=T)


###Here the script in case you prefer (or your reviewer prefers) nmds
# library("RColorBrewer")
# NMDS<-metaMDS(raOTU, distance="bray")
# 
# col8=brewer.pal(12, "Paired")
# colNew<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
# 
# ordiplot(NMDS,type="n")

# orditorp(NMDS,display="sites",cex=0.8,air=0.01, col=colNew[variables$matrix])


MP<-subset(raOTU, variables$matrix=="MP")  #Subset the dataset to only get the MP samples
VMP<-subset(variables, variables$matrix=="MP") #Subset the variables too
water<-subset(raOTU, variables$matrix=="water") #Subset water 
Vwater<-subset(variables, variables$matrix=="water")

betaMP<-vegdist(MP, method="bray") #make betadiversity only for MP & water
adonisMP<-adonis(betaMP~VMP$desin*VMP$date)

betaw<-vegdist(water, method="bray")
adonisW<-adonis(betaw~Vwater$desin*Vwater$date)

adonisMP
adonisW

#Compare the distance matrix of the miroplastic 
rownames(MP)<-VMP$place_time_rep  #set the same names for the MP and water samples from the same place and time (corresponding samples)
rownames(water)<-Vwater$place_time_rep
betaMP2<-vegdist(MP, method="bray")
betaw2<-vegdist(water, method="bray")
par(mfrow=c(2,1))
plot(hclust(betaMP2), main="bray curtis Microplastic", hang=-1)
plot(hclust(betaw2), main="bray curtis water", hang=-1)
mantel<-mantel(betaMP2,betaw2) #Test whether the two distance matrixes are correlated
mantel

write("mantel statistics: MP vs water, r and p value", "stats.csv", append=T)
write(mantel$statistic, "stats.csv", append=T)
write(mantel$signif, "stats.csv", append=T)


####TAXONOMY
library("reshape2")

#taxa as expample on family level only
aggT<-do.call("rbind", as.list(by(traOTU[,], taxa$f, colSums))) #make the sum of the OTU table based on the taxa$f column meaning the family
alldataT<-as.data.frame(rowSums(aggT)) #calculated the total abundace of each family

rareST4<-subset(aggT, alldataT<2000) #make a dataframe with less abundant families
abundST4<-subset(aggT, alldataT>=2000) #...and abundant families
sumrareST4<-as.data.frame(colSums(rareST4)) # sum of rare familie per sample
sumrareST4<-t(sumrareST4)
row.names(sumrareST4)<-"other"  #name them other
colnames(sumrareST4)<-colnames(abundST4)
allST4<-as.data.frame(rbind(abundST4, sumrareST4)) 
#Plot as heat map
mST<-as.matrix(allST4) #trasfer to matrix
heatmap(mST, col=blues9)


#INDVAL ON GENUS
# ###INDVAL
library("indicspecies")
taggT<-as.data.frame(t(aggT))
indval<-multipatt(taggT, variables$matrix) #find inicator species (here families) for the two treatments
summary(indval, alpha=0.005) #summary but show only with p >=0.005

outInd<-capture.output(summary(indval, indvalcomp=TRUE))
write.csv(outInd, "indval.csv")

indPlot<-indval$sign[(indval$sign$p.value<=0.005)&(!is.na(indval$sign$p.value)),] #subset dataframe families where the p value is below 0.005 but not NA


aggIND<-as.data.frame(aggT[row.names(aggT)%in%row.names(indPlot),]) #subset the abundace dataframe to the significant families
aggINDabu<-aggIND[rowSums(aggIND)>=100,] #Only use families with more than 100 reads

taggIND<-as.data.frame(t(aggIND)) #transpose
plot(taggIND[,3]~variables$matrix, ylab="# reads", xlab="matrix", main="Campylobacteraceae") #plot one family

#Make a loop for all the important families: make and export a plot for each of them based on their position in the dataframe
for (i in 1:15) {
  pdf( file= paste("C:/Users/Ester/Desktop/CNR/WWTP/Microplastic/",as.character(i),"FAM_Invdal.pdf",sep=''),height=5,width=2.5)
  plot(taggIND[,i]~variables$matrix, xlab=NULL, ylab="# reads", main=colnames(taggIND)[i])
  graphics.off()
}


###Random Forest analysis###
#Machine learning: How predictable are the communities and which families are the most important in defining them?

library("randomForest")

group<-variables$matrix #How well can the communities can be learned?
predictorMPF<-randomForest(y=as.factor(group), x=taggT, ntree=1000, importance=TRUE) #Try to predict the communities
predictorMPF
pre<-capture.output(predictorMPF)
write("random forest prediction:", "stats.csv", append=T)
write(pre, "stats.csv", append=T)

varImpPlot(predictorMPF)
predictorMPF$importance

group<-variables$desin
predictorD<-randomForest(y=as.factor(group), x=taggT, ntree=1000, importance=TRUE) #
predictorD

write.csv(predictor$importance, "meandecrease.csv")


   


