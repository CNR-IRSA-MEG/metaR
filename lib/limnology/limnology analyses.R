###########################################################################

## The use of the term ‘limnology’ and its scientometrics consequences ####
## Diego Fontaneto, Alejandro Martínez, Stefano Mammola, Aldo Marchetto ###

###########################################################################

##### Statistical analyses ################################################


# Working directory -------------------------------------------------------

setwd("") # <- change me :)

# Loading R packages ------------------------------------------------------

# library("car")
library("emmeans")
library("gam")
library("ggplot2")
# library("MASS")
# library('nlme')
library("performance")
library("plyr")
library("RColorBrewer")

# Load the dataset --------------------------------------------------------

db2 <- read.csv(file="Data_limnology.csv", dec='.', header=TRUE, as.is=T)

names(db2)
summary(db2)

# replace names of 'keyword.search'

db2$keyword.search <- plyr::revalue(db2$keyword.search,
								c("lake_ecolog"="lake ecology",
								"limnolog"="limnology",
								"hydrobiolog"="hydrobiology",
								"oceanograph"="oceanography"))

##### ANALYSES ############################################################


# number of authors -------------------------------------------------------

db2_aut <- db2[db2$N_aut>0,]
db2_aut$keyword.search <- as.factor(db2_aut$keyword.search)

m_authors <- gam::gam(N_aut ~ keyword.search*s(Publication.Year),
			 family=poisson, data=db2_aut)

performance::check_model(m_authors)

# m_authors <- MASS::glm.nb(N_aut ~ keyword.search*Publication.Year,
# 			 data=db2_aut)
# car::Anova(m_authors)

summary(m_authors)

pairs(emmeans::emmeans(m_authors, ~keyword.search*s(Publication.Year)),
	simple="keyword.search")

ggplot2::ggplot(db2_aut, aes(x=Publication.Year, y=N_aut)) + 
	theme_light() +
	geom_point(colour='grey') + 
	geom_smooth(aes(colour=keyword.search)) +
	scale_y_continuous(trans='log2') +
	scale_x_continuous(breaks = c(seq(from=1960, to=2020, by=10)), 
         labels=c("1960","1970","1980","1990","2000","2010","2020")) +
    labs(x=NULL, y="Number of authors", colour="topic") +
    scale_color_brewer(palette="Paired")


# number of references ----------------------------------------------------

db2_ref <- db2[db2$Cited.Reference.Count>0,]
db2_ref$keyword.search <- as.factor(db2_ref$keyword.search)

# m_refs <- MASS::glm.nb(Cited.Reference.Count ~ 
#		keyword.search*Publication.Year, data=db2_ref)
# car::Anova(m_refs)

m_refs <- gam::gam(Cited.Reference.Count ~ 
			keyword.search*s(Publication.Year),
			family=poisson, data=db2_aut)
performance::check_model(m_refs)

summary(m_refs)

pairs(emmeans::emmeans(m_refs, ~keyword.search*s(Publication.Year)),
	simple="keyword.search")

ggplot2::ggplot(db2_ref, aes(x=Publication.Year, 
	y=Cited.Reference.Count)) + 
	theme_light() +
	geom_point(colour='grey') + 
	geom_smooth(aes(colour=keyword.search)) +
	scale_y_continuous(trans='log2') +
	scale_x_continuous(breaks = c(seq(from=1960, to=2020, by=10)), 
         labels=c("1960","1970","1980","1990","2000","2010","2020")) +
    labs(x=NULL, y="Number of references") +
    scale_color_brewer(palette="Paired")

db2_ref90 <- db2_ref[db2_ref$Publication.Year>1989,]
db2_ref90$keyword.search <- as.factor(db2_ref90$keyword.search)
m_refs90 <- gam::gam(Cited.Reference.Count ~ 
		keyword.search*s(Publication.Year), family=poisson, data=db2_ref90)
performance::check_model(m_refs90)
summary(m_refs90)
pairs(emmeans::emmeans(m_refs90, ~keyword.search*s(Publication.Year)),
	simple="keyword.search")

# number of pages ---------------------------------------------------------

db2_pag <- db2[db2$Number.of.Pages>0,]
db2_pag$keyword.search <- as.factor(db2_pag$keyword.search)

# m_pages <- MASS::glm.nb(Number.of.Pages ~
#		 keyword.search*Publication.Year, data=db2_pag)
# car::Anova(m_pages)

m_pages <- gam::gam(Number.of.Pages ~ keyword.search*s(Publication.Year),
			 family=poisson, data=db2_pag)

performance::check_model(m_pages)

summary(m_pages)

pairs(emmeans::emmeans(m_pages, ~keyword.search*s(Publication.Year)),
	simple="keyword.search")

ggplot2::ggplot(db2_pag, aes(x=Publication.Year, y=Number.of.Pages)) + 
	theme_light() +
	geom_point(colour='grey') + 
	geom_smooth(aes(colour=keyword.search)) +
	scale_y_continuous(trans='log2') +
	scale_x_continuous(breaks = c(seq(from=1960, to=2020, by=10)), 
         labels=c("1960","1970","1980","1990","2000","2010","2020")) +
    labs(x=NULL, y="Number of pages") +
    scale_color_brewer(palette="Paired")


db2_pag90 <- db2_pag[db2_pag$Publication.Year>1989,]
db2_pag90$keyword.search <- as.factor(db2_pag90$keyword.search)
m_pages90 <- gam::gam(Number.of.Pages ~ keyword.search*s(Publication.Year),
			 family=poisson, data=db2_pag90)
performance::check_model(m_pages90)
summary(m_pages90)
pairs(emmeans::emmeans(m_pages90, ~keyword.search*Publication.Year),
	simple="keyword.search")


# number of citations -----------------------------------------------------

db2_cit <- db2[db2$Times.Cited..All.Databases>=0,]
db2_cit$keyword.search <- as.factor(db2_cit$keyword.search)

# m_citations <- MASS::glm.nb(Times.Cited..All.Databases ~ 
#		keyword.search*Publication.Year, data=db2_cit)
# car::Anova(m_citations)

m_citations <- gam::gam(Times.Cited..All.Databases ~ 
		keyword.search*s(Publication.Year), family=poisson, data=db2_cit)

performance::check_model(m_citations)

summary(m_citations)

pairs(emmeans::emmeans(m_citations, ~keyword.search*s(Publication.Year)),
	simple="keyword.search")

ggplot2::ggplot(db2_cit, aes(x=Publication.Year, 
		y=Times.Cited..All.Databases)) + 
	theme_light() +
	geom_point(colour='grey') + 
	geom_smooth(aes(colour=keyword.search)) +
	scale_y_continuous(trans='log2') +
	scale_x_continuous(breaks = c(seq(from=1960, to=2020, by=10)), 
         labels=c("1960","1970","1980","1990","2000","2010","2020")) +
    labs(x=NULL, y="Number of citations", colour="topic") +
    scale_color_brewer(palette="Paired")


# readability abstract ---------------------------------------------------

boxplot(db2$FRE_abs)
db2_rea <- db2[db2$FRE_abs > -80,]
length(db2$FRE_abs)-length(db2_rea$FRE_abs)
db2_rea$keyword.search <- as.factor(db2_rea$keyword.search)

# m_reada <- lm(FRE_abs ~ keyword.search*Publication.Year,
#			 data=db2_rea)
# car::Anova(m_reada)

m_reada <- gam::gam(FRE_abs ~ keyword.search*s(Publication.Year), 	
		data=db2_aut)

performance::check_model(m_reada)

summary(m_reada)

pairs(emmeans::emmeans(m_reada, ~keyword.search*s(Publication.Year)),
	simple="keyword.search")

ggplot2::ggplot(db2_rea, aes(x=Publication.Year, y=FRE_abs)) + 
	theme_light() +
	geom_point(colour='grey') + 
	geom_smooth(aes(colour=keyword.search)) +
	scale_x_continuous(breaks = c(seq(from=1960, to=2020, by=10)), 
         labels=c("1960","1970","1980","1990","2000","2010","2020")) +
    labs(x=NULL, y="Readability of abstract", colour="topic") +
    scale_color_brewer(palette="Paired")

db2_rea90 <- db2_rea[db2_rea$Publication.Year>1989,]
db2_rea90$keyword.search <- as.factor(db2_rea90$keyword.search)
m_reada90 <- gam::gam(FRE_abs ~ keyword.search*s(Publication.Year),
			 data=db2_rea90)
performance::check_model(m_reada90)
summary(m_reada90)
pairs(emmeans::emmeans(m_reada90, ~keyword.search*s(Publication.Year)),
	simple="keyword.search")
	
	
# readability title ------------------------------------------------------

boxplot(db2$FRE_tit)
db2_ret <- db2[db2$FRE_tit > -200,]
length(db2$FRE_tit)-length(db2_ret$FRE_tit)
db2_ret$keyword.search <- as.factor(db2_ret$keyword.search)

# m_readt <- lm(FRE_tit ~ keyword.search*Publication.Year,
#			 data=db2_ret)
# car::Anova(m_readt)
			 
m_readt <- gam::gam(FRE_tit ~ keyword.search*s(Publication.Year), 	
		data=db2_aut)

performance::check_model(m_readt)

summary(m_readt)

pairs(emmeans::emmeans(m_readt, ~keyword.search*s(Publication.Year)),
	simple="keyword.search")

ggplot2::ggplot(db2_ret, aes(x=Publication.Year, y=FRE_tit)) + 
	theme_light() +
	geom_point(colour='grey') + 
	geom_smooth(aes(colour=keyword.search)) +
	scale_x_continuous(breaks = c(seq(from=1960, to=2020, by=10)), 
         labels=c("1960","1970","1980","1990","2000","2010","2020")) +
    labs(x=NULL, y="Readability of title", colour="topic") +
    scale_color_brewer(palette="Paired")

db2_ret90 <- db2_ret[db2_ret$Publication.Year>1989,]
db2_ret90$keyword.search <- as.factor(db2_ret90$keyword.search)
m_readt90 <- gam::gam(FRE_tit ~ keyword.search*s(Publication.Year),
			 data=db2_ret90)
performance::check_model(m_readt90)
summary(m_readt90)
pairs(emmeans::emmeans(m_readt90, ~keyword.search*s(Publication.Year)),
	simple="keyword.search")


# citation trends --------------------------------------------------------

library("mgcv") # to load the mgcv library for GAMM
library("itsadug") # to load the mgcv library for Wald test

db2_cit$age <- 2021-db2_cit$Publication.Year # to calculate age of papers

m_cits_res <- mgcv::gam(Times.Cited..All.Databases ~ 
	s(age), family=poisson, data=db2_cit) # to model the effect of age

plot(m_cits_res, se=T)

db2_cit$citation_res <- resid(m_cits_res, type = "pearson") # to extract the residual of citations

db2_cit$cit_res_log <- log((db2_cit$citation_res-min(db2_cit$citation_res))+1) # logarithmic transformation

db2_cit$journal <- as.factor(db2_cit$Source.Title)
db2_cit$keyword.search <- as.factor(db2_cit$keyword.search)

db2_cit_use <- db2_cit[!is.na(db2_cit$FRE_abs),] # remove rows with NA in FRE_abs

dim(db2_cit)
dim(db2_cit_use)

# 1. with gamm
model_citations <- mgcv::gamm(cit_res_log ~ 
		keyword.search +
		s(N_aut, by=keyword.search) + 
		s(FRE_abs, by=keyword.search) + 
		s(Cited.Reference.Count, by=keyword.search) + 
		s(Number.of.Pages, by=keyword.search), 
		random=list(journal=~1), 
		data=db2_cit_use)

summary(model_citations$gam)
anova(model_citations$gam)

wald_gam(model_citations$gam)

boxplot(cit_res_log~keyword.search, 
		data= db2_cit_use,
		xlab="Topic",
		ylab="(Residual) number of citations")

ggplot2::ggplot(db2_cit_use, aes(x=N_aut, y=cit_res_log)) +
	theme_light() +
	geom_smooth(aes(colour=keyword.search)) +
    labs(x="Number of authors", y="(Residual) number of citations",
    	colour="topic") +
    scale_color_brewer(palette="Paired")

ggplot2::ggplot(db2_cit_use, aes(x=FRE_abs, y=cit_res_log)) +
	theme_light() +
	geom_smooth(aes(colour=keyword.search)) +
    labs(x="Abstract readability", y="(Residual) number of citations",
    	colour="topic") +
    scale_color_brewer(palette="Paired")
  
ggplot2::ggplot(db2_cit_use, 
		aes(x=Cited.Reference.Count, y=cit_res_log)) +
	theme_light() +
	geom_smooth(aes(colour=keyword.search)) +
    labs(x="Number of cited references", 
    	y="(Residual) number of citations",
    	colour="topic") +
    scale_color_brewer(palette="Paired")

ggplot2::ggplot(db2_cit_use, aes(x=Number.of.Pages, y=cit_res_log)) +
	theme_light() +
	geom_smooth(aes(colour=keyword.search)) +
    labs(x="Number of pages", y="(Residual) number of citations",
    	colour="topic") +
    scale_color_brewer(palette="Paired")


# 2. with lme

model_citations_LME <- nlme::lme(cit_res_log ~ 
		N_aut*keyword.search + 
		FRE_abs*keyword.search + 
		Cited.Reference.Count*keyword.search + 
		Number.of.Pages*keyword.search, 
		random=list(journal=~1), 
		data=db2_cit_use, method='REML')

car::Anova(model_citations_LME)

model_citations_LME2 <- nlme::lme(cit_res_log ~ 
		N_aut + keyword.search +
		FRE_abs*keyword.search + 
		Cited.Reference.Count + 
		Number.of.Pages*keyword.search, 
		random=list(journal=~1), 
		data=db2_cit_use, method='REML')

car::Anova(model_citations_LME2)



