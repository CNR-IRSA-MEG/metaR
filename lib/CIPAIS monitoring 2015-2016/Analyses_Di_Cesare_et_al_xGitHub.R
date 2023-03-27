## ------------------------------------------------------------------------
## 'ADD TITLE' 
## ------------------------------------------------------------------------

# Di Cesare, A. et al.

# R script to generate the analysis

# Loading R packages ------------------------------------------------------

library("dplyr")
library("glmmTMB")
library("parameters")
library("performance")
library("ggfortify")
library("ggcorrplot")
library("tidylog")
library("tidyverse") 

# Setting working directory -----------------------------------------------

setwd("/Users/your_path/your_path/your_path/your_path") #Change

# Loading the database  ---------------------------------------------------

db <- read.csv2(file = "db.csv", sep = '\t', dec = ',', header = TRUE, as.is = FALSE)

str(db)

# Re-organizing the data to get correlations ------------------------------

db$Date_num <- as.factor(db$Date_num)

db3 <- db[db$Place == "LM Ghiffa", ] ; db3 = droplevels(db3) #all sites except the Lake
db2 <- db[db$Place != "LM Ghiffa", ] ; db2 = droplevels(db2) #only the lake

db_ghiffa <- data.frame(Date = rep(db3$Date_num,5),
                        Type2 = c(rep("sul2.16S",nrow(db3)),
                                rep("tetA.16S",nrow(db3)), 
                                rep("ermB.16S",nrow(db3)),
                                rep("qnrS.16S",nrow(db3)),
                                rep("blaCTXM.16S",nrow(db3))
                                ),
                        Val2 = c(db3$sul2.16S,
                                 db3$tetA.16S,
                                 db3$ermB.16S,
                                 db3$qnrS.16S,
                                 db3$blaCTXM.16S
                                 ))

#creating the variable for left_join
db_ghiffa$DateType <- paste(db_ghiffa$Type,db_ghiffa$Date,sep="_")

#Assembling the database
for(i in 1:nlevels(db2$Place)){ 
  
  db_i <- db2[db2$Place == levels(db2$Place)[i],]
  
  db_place = data.frame(Date = rep(db_i$Date_num,5),
                        Place = rep(levels(db2$Place)[i],nrow(db_i)*5),
                        Type = c(rep("sul2.16S",nrow(db_i)),
                                  rep("tetA.16S",nrow(db_i)), 
                                  rep("ermB.16S",nrow(db_i)),
                                  rep("qnrS.16S",nrow(db_i)),
                                  rep("blaCTXM.16S",nrow(db_i))),
                        Val = c(db_i$sul2.16S,
                                 db_i$tetA.16S,
                                 db_i$ermB.16S,
                                 db_i$qnrS.16S,
                                 db_i$blaCTXM.16S
                                 ))
  
  db_place$DateType <- paste(db_place$Type,db_place$Date,sep="_")
  
  db_place <- db_place  %>%  left_join(db_ghiffa, by = "DateType")
  
  if(i > 1)
    db_cor <- rbind(db_cor,db_place)
  else
      db_cor <- db_place 

  }

str(db_cor)

db_cor <- db_cor %>% dplyr::mutate_if(is.character, as.factor)

db_cor <- db_cor %>% dplyr::select(-c("DateType","Date.y","Type2"))
db_cor$Date.x <- as.numeric(db_cor$Date.x)

# -------------------------------------------------------------------------
# Correlation -------------------------------------------------------------
# -------------------------------------------------------------------------

#Calculating correlations
for(i in 1:nlevels(db_cor$Place)){ 
  
  db_i <- db_cor[db_cor$Place == levels(db_cor$Place)[i],] ; db_i = droplevels(db_i)
  
  result_i <- by(db_i[,c("Val","Val2")], db_i$Type, function(x) {cor(x$Val, x$Val2)})
  
  result_i <- data.frame(as.matrix(result_i)) %>% rownames_to_column("Type")
  colnames(result_i)[2] = "Correlation"
  
  result_i <- data.frame(result_i, Site = rep(levels(db_i$Place),nrow(result_i)))
  
  if(i > 1)
    result <- rbind(result,result_i)
  else
    result <- result_i 
}
warnings()

#Store in a table
correlation_table <- na.omit(result) %>% arrange(Type) ; rm(result_i,results,i)

#Save
write.csv(correlation_table,"correlation_table.csv")

# -------------------------------------------------------------------------
# Relationships -----------------------------------------------------------
# -------------------------------------------------------------------------

#Data exploration
db$Date_num <- as.numeric(db$Date_num)

colnames(db)
#collinearity

db$X.P.reac.mg.L   <- as.numeric(db$X.P.reac.mg.L)
db$X.P..total.mg.L <- as.numeric(db$X.P..total.mg.L)

db %>% dplyr::select(cell.mL,
                     aggregates.mL,
                     microcolonies.mL,
                     X.N.NO3mg.L,
                     X.N.NH4.mg.L,
                     X.P.reac.mg.L,
                     X.P..total.mg.L,
                     total.N.mg.L,
                     TOC.mg.L) %>% 
  GGally::ggpairs() + theme_classic()

db %>% dplyr::select(cell.mL,
                     total.N.mg.L,
                     X.N.NH4.mg.L,
                     Date_num) %>% 
  GGally::ggpairs() + theme_classic()

# Assembling a new dataset with uncollinear variables
db_glm <- db %>%  dplyr::select(cell.mL,
                                total.N.mg.L,
                                X.N.NH4.mg.L,
                                X.P..total.mg.L,
                                Date_num,
                                Place) %>% 
  mutate_if(is.numeric, scale)

db_glm <- data.frame(db %>% dplyr::select(ends_with(".P.A"),
                                          ends_with(".16S")), db_glm)

#Balance of Presence/absence data
db_glm$sul2.P.A %>% table #highly unbalanced!
db_glm$tetA.P.A %>% table #highly unbalanced!
db_glm$ermB.P.A %>% table
db_glm$qnrS.P.A %>% table
db_glm$blaCTXM.P.A %>% table

# Distirbution of concentrations

ggplot(db_glm, aes(x = sul2.16S)) + 
  geom_dotplot(binaxis='x', stackdir='center', dotsize = 0.5) + 
  theme_classic()

ggplot(db_glm, aes(x = tetA.16S)) + 
  geom_dotplot(binaxis='x', stackdir='center', dotsize = 0.5) + 
  theme_classic()

ggplot(db_glm, aes(x = ermB.16S)) + 
  geom_dotplot(binaxis='x', stackdir='center', dotsize = 0.5) + 
  theme_classic()

ggplot(db_glm, aes(x = qnrS.16S)) + 
  geom_dotplot(binaxis='x', stackdir='center', dotsize = 0.5) + 
  theme_classic()

#Outliers in all of them!

# Presence absence models -------------------------------------------------

formula <- "Date_num + cell.mL + total.N.mg.L + X.N.NH4.mg.L + X.P..total.mg.L + (1 | Place)"

formula.sul2.P.A <- as.formula(paste0("sul2.P.A ~ ",formula)) #unbalanced!
formula.tetA.P.A <- as.formula(paste0("tetA.P.A ~ ",formula)) #unbalanced!
formula.ermB.P.A <- as.formula(paste0("ermB.P.A ~ ",formula))
formula.qnrS.P.A <- as.formula(paste0("qnrS.P.A ~ ",formula))
formula.blaCTXM.P.A <- as.formula(paste0("blaCTXM.P.A ~ ",formula))

M1.P.A <- glmmTMB::glmmTMB(formula.sul2.P.A,
                       family = binomial(link = "cloglog"),
                       data = db_glm)

M2.P.A <- glmmTMB::glmmTMB(formula.tetA.P.A,
                       family = binomial(link = "cloglog"),
                       data = db_glm)

M3.P.A <- glmmTMB::glmmTMB(formula.ermB.P.A,
                       family = binomial(link = "cloglog"), 
                       data = db_glm)

M4.P.A <- glmmTMB::glmmTMB(formula.qnrS.P.A,
                       family = binomial(link = "cloglog"), 
                       data = db_glm)

M5.P.A <- glmmTMB::glmmTMB(formula.blaCTXM.P.A,
                           family = binomial(link = "cloglog"), 
                           data = db_glm)

#Model validation
performance::check_model(M1.P.A) #so so (consequence of low sample size and unbalnce of 0/1)
performance::check_model(M2.P.A) #so so (consequence of low sample size and unbalnce of 0/1)
performance::check_model(M3.P.A) #so so (consequence of low sample size and unbalnce of 0/1)
performance::check_model(M4.P.A) #so so (consequence of low sample size and unbalnce of 0/1)
performance::check_model(M5.P.A) #so so (consequence of low sample size and unbalnce of 0/1)

#Summary tables
(par.M1 <- parameters::parameters(M1.P.A))
(par.M2 <- parameters::parameters(M2.P.A))
(par.M3 <- parameters::parameters(M3.P.A))
(par.M4 <- parameters::parameters(M4.P.A))
(par.M5 <- parameters::parameters(M5.P.A))

table.M1 <- par.M1 %>% dplyr::select(Parameter,
                                     Effects,
                                     Beta = Coefficient,
                                     SE,
                                     CI_low,
                                     CI_high,
                                     z,
                                     p) %>% 
  data.frame() %>% 
  mutate_if(is.numeric, ~ round(.,3)) ; rm(par.M1)

table.M1 <- table.M1[table.M1$Effects == "fixed",] %>% 
  dplyr::select(-c(Effects)) %>% 
  na.omit()

table.M2 <- par.M2 %>% dplyr::select(Parameter,
                                     Effects,
                                     Beta = Coefficient,
                                     SE,
                                     CI_low,
                                     CI_high,
                                     z,
                                     p) %>% 
  data.frame() %>% 
  mutate_if(is.numeric, ~ round(.,3)) ; rm(par.M2)

table.M2 <- table.M2[table.M2$Effects == "fixed",] %>% 
  dplyr::select(-c(Effects)) %>% 
  na.omit()

table.M3 <- par.M3 %>% dplyr::select(Parameter,
                                     Effects,
                                     Beta = Coefficient,
                                     SE,
                                     CI_low,
                                     CI_high,
                                     z,
                                     p) %>% 
  data.frame() %>% 
  mutate_if(is.numeric, ~ round(.,3)) ; rm(par.M3)

table.M3 <- table.M3[table.M3$Effects == "fixed",] %>% 
  dplyr::select(-c(Effects)) %>% 
  na.omit()

table.M4 <- par.M4 %>% dplyr::select(Parameter,
                                     Effects,
                                     Beta = Coefficient,
                                     SE,
                                     CI_low,
                                     CI_high,
                                     z,
                                     p) %>% 
  data.frame() %>% 
  mutate_if(is.numeric, ~ round(.,3)) ; rm(par.M4)

table.M4 <- table.M4[table.M4$Effects == "fixed",] %>% 
  dplyr::select(-c(Effects)) %>% 
  na.omit()

table.M5 <- par.M5 %>% dplyr::select(Parameter,
                                     Effects,
                                     Beta = Coefficient,
                                     SE,
                                     CI_low,
                                     CI_high,
                                     z,
                                     p) %>% 
  data.frame() %>% 
  mutate_if(is.numeric, ~ round(.,3)) ; rm(par.M5)

table.M5 <- table.M5[table.M5$Effects == "fixed",] %>% 
  dplyr::select(-c(Effects)) %>% 
  na.omit()

table.M <- cbind(Model = c(rep("sul2 [Presence/Absence]",nrow(table.M1)),
                           rep("tetA [Presence/Absence]",nrow(table.M2)),
                           rep("ermB [Presence/Absence]",nrow(table.M3)),
                           rep("qnrS [Presence/Absence]",nrow(table.M4)),
                           rep("blaCTXM [Presence/Absence]",nrow(table.M5))),
                 rbind(table.M1,table.M2,table.M3,table.M4,table.M5)) ; rm(table.M1,table.M2,table.M3,table.M4)

table.M$Parameter <- as.factor(as.character(table.M$Parameter))

#rename
levels(table.M$Parameter) <- c("Intercept", 
                               "total cell number", 
                               "Sampling date", 
                               "total N (TN)", 
                               "N-NH4",
                               "total P (TP)")

table.M$Parameter <- factor(table.M$Parameter, c("Intercept", 
                                                 "total cell number", 
                                                 "Sampling date", 
                                                 "total N (TN)", 
                                                 "N-NH4",
                                                 "total P (TP)")) #Sort

# Save
write.csv(table.M,"results_P_A_models.csv")

# Plot
table.plot.PA <- table.M[table.M$Parameter != "Intercept",] ; table.plot.PA = droplevels(table.plot.PA)

sign.PA <- ifelse(table.plot.PA$p > 0.05, "", ifelse(table.plot.PA$p > 0.01," *", " **")) #Significance

(PA.forest_plot <- 
    table.plot.PA %>%
    ggplot2::ggplot(aes(x = Beta, y = Parameter)) + 
    facet_wrap(. ~ Model, nrow = 3, ncol = 2, scales='free') +  
    geom_vline(lty = 3, size = 0.5, col = "grey50", xintercept = 0) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0, col = "grey15")+
    geom_point(size = 2, pch = 21, col = "grey15", fill = "grey30") +
    geom_text(label = paste0(round(table.plot.PA$Beta, 3), sign.PA, sep = "  "), 
              vjust = - 1, size = 3) +
    labs(x = expression(paste("Estimated beta" %+-% "95% Confidence interval")),
         y = NULL) +
    
    theme_classic() + 
    theme(axis.line = element_line(),
                            strip.background = element_blank(),
                            strip.text = element_text(face = "bold"),
                            panel.margin = unit(.9, "lines")))

# Continuous models -------------------------------------------------

formula <- "Date_num + cell.mL + total.N.mg.L + X.N.NH4.mg.L +  X.P..total.mg.L + (1 | Place)"

#too many zeros!
formula.sul2.16S <- as.formula(paste0("sul2.16S ~ ",formula))
formula.tetA.16S <- as.formula(paste0("tetA.16S ~ ",formula))
formula.ermB.16S <- as.formula(paste0("ermB.16S ~ ",formula))
formula.qnrS.16S <- as.formula(paste0("qnrS.16S ~ ",formula))

M1.16S <- glmmTMB::glmmTMB(formula.sul2.16S,
                           family = gaussian, 
                           data = db_glm[db_glm$sul2.16S<0.02,]) #removing outlier

M2.16S <- glmmTMB::glmmTMB(formula.tetA.16S,
                           family = gaussian, 
                           data = db_glm[db_glm$tetA.16S<0.02,])

M3.16S <- glmmTMB::glmmTMB(formula.ermB.16S,
                           family = gaussian, 
                           data = db_glm[db_glm$ermB.16S<0.02,])

M4.16S <- glmmTMB::glmmTMB(formula.qnrS.16S,
                           family = gaussian, 
                           data = db_glm[db_glm$qnrS.16S<0.005,])

#Model validation
performance::check_model(M1.16S) #ok
performance::check_model(M2.16S) #so so
performance::check_model(M3.16S) #so so
performance::check_model(M4.16S) #so so

#Summary tables
(par.M1 <- parameters::parameters(M1.16S))
(par.M2 <- parameters::parameters(M2.16S))
(par.M3 <- parameters::parameters(M3.16S))
(par.M4 <- parameters::parameters(M4.16S))

table.M1 <- par.M1 %>% dplyr::select(Parameter,
                                     Effects,
                                     Beta = Coefficient,
                                     SE,
                                     CI_low,
                                     CI_high,
                                     z,
                                     p) %>% 
  data.frame() %>% 
  mutate_if(is.numeric, ~ round(.,6)) ; rm(par.M1)

table.M1 <- table.M1[table.M1$Effects == "fixed",] %>% 
  dplyr::select(-c(Effects)) %>% 
  na.omit()

table.M2 <- par.M2 %>% dplyr::select(Parameter,
                                     Effects,
                                     Beta = Coefficient,
                                     SE,
                                     CI_low,
                                     CI_high,
                                     z,
                                     p) %>% 
  data.frame() %>% 
  mutate_if(is.numeric, ~ round(.,6)) ; rm(par.M2)

table.M2 <- table.M2[table.M2$Effects == "fixed",] %>% 
  dplyr::select(-c(Effects)) %>% 
  na.omit()

table.M3 <- par.M3 %>% dplyr::select(Parameter,
                                     Effects,
                                     Beta = Coefficient,
                                     SE,
                                     CI_low,
                                     CI_high,
                                     z,
                                     p) %>% 
  data.frame() %>% 
  mutate_if(is.numeric, ~ round(.,6)) ; rm(par.M3)

table.M3 <- table.M3[table.M3$Effects == "fixed",] %>% 
  dplyr::select(-c(Effects)) %>% 
  na.omit()

table.M4 <- par.M4 %>% dplyr::select(Parameter,
                                     Effects,
                                     Beta = Coefficient,
                                     SE,
                                     CI_low,
                                     CI_high,
                                     z,
                                     p) %>% 
  data.frame() %>% 
  mutate_if(is.numeric, ~ round(.,6)) ; rm(par.M4)

table.M4 <- table.M4[table.M4$Effects == "fixed",] %>% 
  dplyr::select(-c(Effects)) %>% 
  na.omit()

table.M <- cbind(Model = c(rep("sul2 [Concentration]",nrow(table.M1)),
                           rep("tetA [Concentration]",nrow(table.M2)),
                           rep("ermB [Concentration]",nrow(table.M3)),
                           rep("qnrS [Concentration]",nrow(table.M4))),
                 rbind(table.M1,table.M2,table.M3,table.M4)) ; rm(table.M1,table.M2,table.M3,table.M4)

table.M$Parameter <- as.factor(as.character(table.M$Parameter))

#rename
levels(table.M$Parameter) <- c("Intercept", 
                               "total cell number", 
                               "Sampling date", 
                               "total N (TN)", 
                               "N-NH4",
                               "total P (TP)")

table.M$Parameter <- factor(table.M$Parameter, c("Intercept", 
                                                 "total cell number", 
                                                 "Sampling date", 
                                                 "total N (TN)", 
                                                 "N-NH4",
                                                 "total P (TP)")) #Sort

# Save
write.csv(table.M,"results_16S_models.csv")

# Plot
table.plot.c <- table.M[table.M$Parameter != "Intercept",] ; table.plot.c = droplevels(table.plot.c)

sign.c <- ifelse(table.plot.c$p > 0.05, "", ifelse(table.plot.c$p > 0.01," *", " **")) #Significance

(c.forest_plot <- 
    table.plot.c  %>%
    ggplot2::ggplot(aes(x = Beta, y = Parameter)) + 
    facet_wrap(. ~ Model, nrow = 2, ncol = 2, scales = "free") +  
    geom_vline(lty = 3, size = 0.5, col = "grey50", xintercept = 0) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0, col = "grey15")+
    geom_point(size = 2, pch = 21, col = "grey15", fill = "grey30") +
    geom_text(label = paste0(round(table.plot.c$Beta, 5), sign.c, sep = "  "), 
              vjust = - 1, size = 3) +
    labs(x = expression(paste("Estimated beta" %+-% "95% Confidence interval")),
         y = NULL) +
    
    theme_classic() + 
    theme(axis.line = element_line(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.margin = unit(.9, "lines")))

# Save Figures ------------------------------------------------------------

#figure S1
pdf(file = "Figure_S1.pdf", width = 12, height = 12)
db_plot = db %>% dplyr::select(cell.mL,
                     aggregates.mL,
                     microcolonies.mL,
                     X.N.NO3mg.L,
                     X.N.NH4.mg.L,
                     total.N.mg.L,
                     TOC.mg.L,
                     X.P.reac.mg.L,
                     X.P..total.mg.L)
colnames(db_plot) = c("total cell number",
                      "number of aggregates", 
                      "number of microcolonies",
                      "N-NO3", "N-NH4", "TN", "TOC", "RP","TP")

db_plot %>%  GGally::ggpairs() + theme_classic()
dev.off()

pdf(file = "Figure_2.pdf", width = 8, height = 7)
Figure_1
dev.off()

pdf(file = "Figure_3.pdf", width = 8, height = 7)
PA.forest_plot
dev.off()

pdf(file = "Figure_4.pdf", width = 8, height = 7)
c.forest_plot
dev.off()
#end

# Correlation with site level variables ----------------------------------

db2 <- read.csv2(file = "var_sito_Gluca_no_dens.csv", sep = ';', dec = ',', header = TRUE, as.is = FALSE)

db_pca = db2[, c(7,9:ncol(db2))]

colnames(db_pca) = c("aggregate", "stream reach", "basin area", "population",
                     "primary industries", "secondary industries", "hotels")

# Principal component analysis
pca <- prcomp(db_pca, scale. = TRUE)

# Plot
(pca.plot <- autoplot(pca, 
                      loadings = TRUE,
                      loadings.colour =  "black",
                      loadings.label = TRUE,
                      loadings.label.colour = "black",
                      data = db2, 
                      fill = 'Site', 
                      size = 4, 
                      shape = 21) + 
    scale_fill_discrete(name=NULL) +
    theme_bw())

db3 <- data.frame(db2[, c(2:7,9:ncol(db2))]) #,PC1 = pca$x[,1], PC2 = pca$x[,2]) 

colnames(db3)[6:ncol(db3)] = c("aggregate", "stream reach", "basin area", "population",
                     "primary industries", "secondary industries", "hotels")

corr <- round(cor(db3), 3)

corr[1:5,1:5] = NA
corr[,6:ncol(corr)] = NA

(cor_plot <- ggcorrplot(corr, hc.order = FALSE,
                        legend.title = "Correlation",
                        type = "lower",
                        outline.col = "white", 
                        insig = "blank"))

pdf(file = "Figure_4.pdf", width = 7, height = 8)
ggpubr::ggarrange(pca.plot,cor_plot,
                  common.legend = FALSE,
                  labels = c("A", "B"),
                  align = "v",
                  ncol=1, nrow=2) 
dev.off()
#end