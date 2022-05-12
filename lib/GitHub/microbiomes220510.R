##### load libraries #####

library(ape)
library(car)
library(ggplot2)
library(vegan)


##### load data #####

roti <- read.csv("rotifer_data_to_use.csv", header=T, as.is=F)
names(roti[,1:15])
summary(roti[,1:10])
head(roti[,1:10])
table(roti$species)

tree <- ape::read.tree('BEAST_combined23.newick')
plot(tree)
tree$tip.label


##### Bray-Curtis distances between 93 samples #####

roti_dist <- vegan::vegdist(roti[11:length(names(roti))],"bray")

# plot the differences

roti_nMDS <- vegan::metaMDS(roti[11:length(names(roti))],"bray")
# plot(roti_nMDS)
data.scores <- as.data.frame(scores(roti_nMDS))
data.scores$species <- roti$species
data.scores$env_cult <- roti$env_cult
data.scores$origin <- roti$origin

nMDS_plot <- ggplot(data.scores, aes(x=NMDS1, y=NMDS2)) + 
    geom_point(size=4, aes(shape=origin, colour=env_cult))+ 
    theme(axis.text.y=element_text(colour="black", size=12, face="bold"), 
    axis.text.x=element_text(colour="black", face="bold", size=12), 
    legend.text=element_text(size=12, face="bold", colour="black"), 
    legend.position="right", axis.title.y=element_text(face="bold", size=14), 
    axis.title.x=element_text(face="bold", size=14, colour="black"), 
    legend.title=element_text(size=14, colour="black", face="bold"), 
    panel.background=element_blank(), panel.border=element_rect(colour="black", fill=NA, size=1.2),
    legend.key=element_blank()) + 
    labs(x="nMDS1", colour="lab/field", y="nMDS2", shape="origin")  + 
    scale_colour_manual(values=c("#ecc6c4","#a74b44")) 

nMDS_plot

# and

plot(hclust(roti_dist, method="average"), labels=roti$species, hang=-1)


##### test for effects on the differences lab/field #####

# one sample for each species

roti_reduced <- droplevels(roti[tapply(1:nrow(roti), roti$species, some, 1),])
dim(roti)
dim(roti_reduced)
roti$species
roti_reduced$species

roti_reduced_dist <- vegan::vegdist(roti_reduced[11:length(names(roti_reduced))],"bray")

effects_reduced <- vegan::adonis(roti_reduced_dist ~ env_cult, data=roti_reduced)
effects_reduced


##### correlation phylogeny-microbiome #####
##### analyses on 23 species, not on 93 samples #####
##### after a random selection of one sample for each species #####

DNA_dist_23 <- ape::cophenetic.phylo(tree)

# check that names of 23 samples are ordered in the same way

# remove spaces and dots from species name in dataset
roti_reduced$species2 <- gsub(" ", "", roti_reduced$species, fixed = TRUE)
roti_reduced$species2 <- gsub(".", "", roti_reduced$species2, fixed = TRUE)
roti_reduced$species2 <- gsub("Synchaetamonopus", "Synchaetasp", roti_reduced$species2, fixed = TRUE)

all(rownames(DNA_dist_23)==roti_reduced$species2)
# which(rownames(DNA_dist_23)!=roti_reduced$species2)

# they are not, thus, reorder and check that orders are fine
roti_ord <- roti_reduced[match(rownames(DNA_dist_23), roti_reduced$species2),]
all(rownames(DNA_dist_23)==roti_ord$species2)

# obtain matrix of Bray-Curtis distances, ordered as in the matrix of phylogenetic distances

roti_ord_dist <- vegan::vegdist(roti_ord[11:(length(names(roti_ord))-1)],"bray")

vegan::mantel(roti_ord_dist,DNA_dist_23)
plot(as.vector(as.dist(DNA_dist_23/max(DNA_dist_23))),as.vector(1-roti_ord_dist), xlab="relative phylogenetic distance", ylab="community similarity")
# plot(vegan::mantel.correlog((1-roti_ord_dist), DNA_dist, cutoff=F, n.class=7))


##### analyses on only environmental samples #####

rotiEnv <- droplevels(subset(roti_ord, env_cult=="environment"))
summary(rotiEnv[,1:10])
dim(roti_ord)
dim(rotiEnv)

rotiEnv_dist <- vegan::vegdist(rotiEnv[11:(length(names(rotiEnv))-1)],"bray")

DNA_distEnv <- subset(DNA_dist_23,roti_ord$env_cult=="environment",roti_ord$env_cult=="environment")
dim(DNA_distEnv)

all(rownames(DNA_distEnv)==rotiEnv$species2)

vegan::mantel(rotiEnv_dist,DNA_distEnv)
plot(as.vector(as.dist(DNA_distEnv/max(DNA_distEnv))),as.vector(1-rotiEnv_dist), xlab="relative phylogenetic distance", ylab="community similarity")

