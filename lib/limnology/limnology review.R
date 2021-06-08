###########################################################################

## The use of the term ‘limnology’ and its scientometrics consequences ####
## Diego Fontaneto, Alejandro Martínez, Stefano Mammola, Aldo Marchetto ###

###########################################################################


## R code to prepare the data:

# Working directory -------------------------------------------------------

setwd("") # <- change me :)

# Loading R packages ------------------------------------------------------

library("bazar")
library("dplyr")
library("tidyr")
library("ggplot2")
library("gridExtra")
library("gtable")
library("grid")
library("jaod")
library("quanteda")

# Functions ---------------------------------------------------------------

## Function for checking if a number is a whole number
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

## Function for calculating the standard error
std <- function(x) sd(x)/sqrt(length(x))

## Function for replacing nan with na
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

## Function for calculating the Flesch reading ease 
FRE_index <- function(wordCount, sentenceCount, syllableCount) {
    FRE <- 206.835 - (1.015 * (wordCount / sentenceCount)) - (84.6 * (syllableCount / wordCount))
  	return(FRE)
}

## Function for  cleaning abstracts
word.cleaner <- function(word.list, remove.punctuation=FALSE, remove.numbers=FALSE, split=FALSE, split.sep= " ") {
  
  #check if right data is provided
  if (class(word.list) != "character"){ 
    word.list <- as.character(word.list)
    warning("The input data as been converted as a character string")
  }
  
  if(split)
    word.list <- strsplit(word.list, split.sep)[[1]]
  
  #all lower case
  word.list <- tolower(as.vector(word.list))

  if(remove.punctuation){
    # remove all punctuation except . , !, and ?
    word.list <- gsub("[^[:alnum:][:space:](?|!|.)]", "\\1", word.list)
    
    # Replace each dot that is in between word characters
    word.list <- gsub("\\b\\.\\b", "\\1", word.list, perl=TRUE)
    
    # Replace each dot that is in between letters
    word.list <- gsub("(?<=\\p{L})\\.(?=\\p{L})", "\\1", word.list, perl=TRUE)
    
    # Replace each dot that is in between word characters, but no in URLs
    word.list <- gsub("(?:ht|f)tps?://\\S*(*SKIP)(*F)|\\b\\.\\b", "\\1", word.list, perl=TRUE)
    
    #remove "et al."
    word.list <- gsub("et al.", "et al", word.list)
    
    #remove "Sic!"
    word.list <- gsub("sic!", "sic", word.list)
    
}
    
  if(remove.numbers)
    word.list <- gsub('[0-9]+', '', word.list)

  #remove extra white spaces
  word.list <- trimws(word.list)
  
  return(word.list)
}

#Example...
test <- "How Are you? Fine, thanks! The paper by Ripple et al.2017 is great. you can download it at www.bioscience.com (Sic!), and enjoy it. v.e.r.y c.o.o.o.l thanks."
word.cleaner(test, remove.punctuation = TRUE, remove.numbers= TRUE, split = FALSE)

## function for counting words, syllabes, and sentences
abstract.counter  <- function(word.list) { 

if (class(word.list) != "character"){ 
  word.list <- as.character(word.list)
  warning("The input data as been converted as a character string")
}
  wordCount <- length(unlist(strsplit(as.character(tokens(word.list, remove_punct = TRUE)), " ")))
  syllableCount <- sum( na.omit(quanteda::nsyllable(tokens(word.list, remove_punct = TRUE))$text1) )
  sentenceCount <- length( unlist(strsplit(word.list, "[.]")) )
  
  output <- list(wordCount,syllableCount,sentenceCount) ; names(output) <- c("wordCount","syllableCount","sentenceCount")
  
  return(output)
  
}

#Example...
warning_insect <- "Here we build on the manifesto 'World Scientists' Warning to Humanity, issued by the Alliance of World Scientists. As a group of conservation biologists deeply concerned about the decline of insect populations, we here review what we know about the drivers of insect extinctions, their consequences, and how extinctions can negatively impact humanity. We are causing insect extinctions by driving habitat loss, degradation, and fragmentation, use of polluting and harmful substances, the spread of invasive species, global climate change, direct overexploitation, and co-extinction of species dependent on other species. With insect extinctions, we lose much more than species. We lose abundance and biomass of insects, diversity across space and time with consequent homogenization, large parts of the tree of life, unique ecological functions and traits, and fundamental parts of extensive networks of biotic interactions. Such losses lead to the decline of key ecosystem services on which humanity depends. From pollination and decomposition, to being resources for new medicines, habitat quality indication and many others, insects provide essential and irreplaceable services. We appeal for urgent action to close key knowledge gaps and curb insect extinctions. An investment in research"

example <- word.cleaner(warning_insect,  remove.punctuation = TRUE, remove.numbers= TRUE, split = FALSE)
example <- abstract.counter(example)

FRE_index(wordCount=example$wordCount, 
          sentenceCount=example$sentenceCount, 
          syllableCount=example$syllableCount) 


# Loading the databases #############################

db <- read.csv(file = "db_all.txt", sep = '\t', dec = '.', header = TRUE)
str(db)
dim(db)

# remove missing years (likely papers still in press in 2021)
dim(db)
db <- db %>% drop_na(Publication.Year)
dim(db)

# remove papers in 2021
dim(db)
db <- db[db$Publication.Year<2021,]
dim(db)

# restrict to journal papers only
dim(db)
db <- db[db$Publication.Type=="J",]
dim(db)


#calculating the number of authors and the FRE readability index for each abstract and title
N_aut <- c()
FRE_abs <- c()
FRE_tit <- c()

#run from here:
for (i in 1:nrow(db)){
  
  db_i <- db[i,]
  
  #Number of authors
  N_aut <- append(N_aut,length(word.cleaner(word.list=as.character(db_i$Authors),
                                            remove.punctuation=FALSE, 
                                            remove.numbers=FALSE, 
                                            split=TRUE, split.sep= ";")) )
  
  #Readability abstract & title
  if(is.na(db_i$Abstract) == TRUE){
    FRE_abs_i <- NA
  }
  
  else {
    clean_abs <- word.cleaner(as.character(db_i$Abstract),  remove.punctuation = TRUE, remove.numbers= TRUE, split = FALSE)
    clean_abs <- abstract.counter(clean_abs)
    FRE_abs_i <- FRE_index(wordCount = clean_abs$wordCount, sentenceCount = clean_abs$sentenceCount, syllableCount = clean_abs$syllableCount)
  }
  
  if(is.na(db_i$Article.Title) == TRUE){
    FRE_tit_i <- NA
  }
  
  else {
    clean_tit <- word.cleaner(as.character(db_i$Article.Title),  remove.punctuation = TRUE, remove.numbers= TRUE, split = FALSE)
    clean_tit <- abstract.counter(clean_tit)
    FRE_tit_i <- FRE_index(wordCount = clean_tit$wordCount, sentenceCount = clean_tit$sentenceCount, syllableCount = clean_tit$syllableCount)
  }
  
  FRE_abs <- append(FRE_abs,FRE_abs_i)
  FRE_tit <- append(FRE_tit,FRE_tit_i)
  
  # Checking the advancing 
  if(is.wholenumber(i/100) == TRUE)
    message(paste("Analyzed", as.character(i), "papers out of", as.character(nrow(db)),sep=" "))

}
#end

# store and save
db2 <- data.frame(db,
				FRE_abs,
				FRE_tit,
				N_aut)
write.csv(db2, "Data_limnology.csv")

###########################################################################
# starting from already prepared data
# db2 <- read.csv(file="Data_limnology.csv", dec='.', header=TRUE, as.is=T)
###########################################################################

# Seeing some trends ------------------------------------------------------

# create datasets for each topic

db_LI <- subset(db2,keyword.search=="limnolog")
db_OC <- subset(db2,keyword.search=="oceanograph")
db_LA <- subset(db2,keyword.search=="lake_ecolog")
db_HY <- subset(db2,keyword.search=="hydrobiolog")

dim(db_LI)
dim(db_OC)
dim(db_LA)
dim(db_HY)


# yearly trends for each topic

db_LI_year <-
  db_LI %>% 
  group_by(Publication.Year) %>% 
  summarise(FRE_mean_abs = mean(FRE_abs, na.rm=T), 
            FRE_se_abs = std(FRE_abs),
            FRE_mean_tit = mean(FRE_tit),
            FRE_se_tit = std(FRE_tit),
            N_aut_mean = mean(N_aut), 
            N_aut_se = std(N_aut),
            N_pages_mean = mean(Number.of.Pages[Number.of.Pages>0]), 
            N_pages_se = std(Number.of.Pages[Number.of.Pages>0]),
            N_citedref_mean = mean(Cited.Reference.Count), 
            N_citedref_se = std(Cited.Reference.Count),
            N_citations_mean = mean(Times.Cited..WoS.Core), 
            N_citations_se = std(Times.Cited..WoS.Core),
            N_papers = length(Publication.Year),
            ) 

db_OC_year <-
  db_OC %>% 
  group_by(Publication.Year) %>% 
  summarise(FRE_mean_abs = mean(FRE_abs, na.rm=T), 
            FRE_se_abs = std(FRE_abs),
            FRE_mean_tit = mean(FRE_tit),
            FRE_se_tit = std(FRE_tit),
            N_aut_mean = mean(N_aut), 
            N_aut_se = std(N_aut),
            N_pages_mean = mean(Number.of.Pages[Number.of.Pages>0]), 
            N_pages_se = std(Number.of.Pages[Number.of.Pages>0]),
            N_citedref_mean = mean(Cited.Reference.Count), 
            N_citedref_se = std(Cited.Reference.Count),
            N_citations_mean = mean(Times.Cited..WoS.Core), 
            N_citations_se = std(Times.Cited..WoS.Core),
            N_papers = length(Publication.Year),
            ) 

db_LA_year <-
  db_LA %>% 
  group_by(Publication.Year) %>% 
  summarise(FRE_mean_abs = mean(FRE_abs, na.rm=T), 
            FRE_se_abs = std(FRE_abs),
            FRE_mean_tit = mean(FRE_tit),
            FRE_se_tit = std(FRE_tit),
            N_aut_mean = mean(N_aut), 
            N_aut_se = std(N_aut),
            N_pages_mean = mean(Number.of.Pages[Number.of.Pages>0]), 
            N_pages_se = std(Number.of.Pages[Number.of.Pages>0]),
            N_citedref_mean = mean(Cited.Reference.Count), 
            N_citedref_se = std(Cited.Reference.Count),
            N_citations_mean = mean(Times.Cited..WoS.Core), 
            N_citations_se = std(Times.Cited..WoS.Core),
            N_papers = length(Publication.Year),
            ) 

db_HY_year <-
  db_HY %>% 
  group_by(Publication.Year) %>% 
  summarise(FRE_mean_abs = mean(FRE_abs, na.rm=T), 
            FRE_se_abs = std(FRE_abs),
            FRE_mean_tit = mean(FRE_tit),
            FRE_se_tit = std(FRE_tit),
            N_aut_mean = mean(N_aut), 
            N_aut_se = std(N_aut),
            N_pages_mean = mean(Number.of.Pages[Number.of.Pages>0]), 
            N_pages_se = std(Number.of.Pages[Number.of.Pages>0]),
            N_citedref_mean = mean(Cited.Reference.Count), 
            N_citedref_se = std(Cited.Reference.Count),
            N_citations_mean = mean(Times.Cited..WoS.Core), 
            N_citations_se = std(Times.Cited..WoS.Core),
            N_papers = length(Publication.Year),
            ) 

dim(db_LI_year)
dim(db_OC_year)
dim(db_LA_year)
dim(db_HY_year)

summary(db_LI_year)
summary(db_OC_year)
summary(db_LA_year)
summary(db_HY_year)


##### PLOTS ###############################################################


# Readability abstract plot -----------------------------------------------

(FRE_a_LI <- ggplot(data=db_LI_year, aes(x=Publication.Year, y=FRE_mean_abs)) + 

   geom_line(aes(x=Publication.Year, y=FRE_mean_abs),linetype = 1,col="blue",alpha=1) + 
   
   #geom_errorbar(aes(ymin=FRE_mean_abs-FRE_se_abs, ymax=FRE_mean_abs+FRE_se_abs), width=.8,col="blue") +
   #geom_point(size = 2,col="blue")+
   
   labs(title= "Limnology",
        x = NULL, 
        y = "Readability [mean]") + #Flesch-Kincaid Readability Ease
   
   scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 
   
   theme_bw()+
   theme(
     legend.position = "none",
     axis.title = element_text(size = 10),
     axis.text.x = element_text(size = 10),
     axis.text.y = element_text(size = 10),
     panel.grid = element_blank(),
     plot.caption = element_text(size = 10, color = "gray50"),
     plot.title = element_text(face="bold", size=12)
   )
)

(FRE_a_OC <- ggplot(data=db_OC_year, aes(x=Publication.Year, y=FRE_mean_abs)) + 

   geom_line(aes(x=Publication.Year, y=FRE_mean_abs),linetype = 1,col="blue",alpha=1) + 
   
   #geom_errorbar(aes(ymin=FRE_mean_abs-FRE_se_abs, ymax=FRE_mean_abs+FRE_se_abs), width=.8,col="blue") +
   #geom_point(size = 2,col="blue")+
   
   labs(title= "Oceanography",
        x = NULL, 
        y = "Readability [mean]") + #Flesch-Kincaid Readability Ease
   
   scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 
   
   theme_bw()+
   theme(
     legend.position = "none",
     axis.title = element_text(size = 10),
     axis.text.x = element_text(size = 10),
     axis.text.y = element_text(size = 10),
     panel.grid = element_blank(),
     plot.caption = element_text(size = 10, color = "gray50"),
     plot.title = element_text(face="bold", size=12)
   )
)

(FRE_a_LA <- ggplot(data=db_LA_year, aes(x=Publication.Year, y=FRE_mean_abs)) + 

   geom_line(aes(x=Publication.Year, y=FRE_mean_abs),linetype = 1,col="blue",alpha=1) + 
   
   #geom_errorbar(aes(ymin=FRE_mean_abs-FRE_se_abs, ymax=FRE_mean_abs+FRE_se_abs), width=.8,col="blue") +
   #geom_point(size = 2,col="blue")+
   
   labs(title= "Lake Ecology",
        x = NULL, 
        y = "Readability [mean]") + #Flesch-Kincaid Readability Ease
   
   scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 
   
   theme_bw()+
   theme(
     legend.position = "none",
     axis.title = element_text(size = 10),
     axis.text.x = element_text(size = 10),
     axis.text.y = element_text(size = 10),
     panel.grid = element_blank(),
     plot.caption = element_text(size = 10, color = "gray50"),
     plot.title = element_text(face="bold", size=12)
   )
)

(FRE_a_HY <- ggplot(data=db_HY_year, aes(x=Publication.Year, y=FRE_mean_abs)) + 

   geom_line(aes(x=Publication.Year, y=FRE_mean_abs),linetype = 1,col="blue",alpha=1) + 
   
   #geom_errorbar(aes(ymin=FRE_mean_abs-FRE_se_abs, ymax=FRE_mean_abs+FRE_se_abs), width=.8,col="blue") +
   #geom_point(size = 2,col="blue")+
   
   labs(title= "Hydrobiology",
        x = NULL, 
        y = "Readability [mean]") + #Flesch-Kincaid Readability Ease
   
   scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 
   
   theme_bw()+
   theme(
     legend.position = "none",
     axis.title = element_text(size = 10),
     axis.text.x = element_text(size = 10),
     axis.text.y = element_text(size = 10),
     panel.grid = element_blank(),
     plot.caption = element_text(size = 10, color = "gray50"),
     plot.title = element_text(face="bold", size=12)
   )
)

# four plots together

FRE_a1 <- ggplotGrob(FRE_a_LI)
FRE_a2 <- ggplotGrob(FRE_a_OC)
FRE_a3 <- ggplotGrob(FRE_a_LA)
FRE_a4 <- ggplotGrob(FRE_a_HY)
grid.arrange(
	FRE_a1, FRE_a2, FRE_a3, FRE_a4,
	nrow=4,
	top="Readability of Abstract"
	)

# Readability title plot --------------------------------------------------

(FRE_t_LI <- ggplot(data=db_LI_year, aes(x=Publication.Year, y=FRE_mean_tit)) + 

   geom_line(aes(x=Publication.Year, y=FRE_mean_tit),linetype = 1,col="blue",alpha=1) + 
   
   geom_errorbar(aes(ymin=FRE_mean_tit-FRE_se_tit, ymax=FRE_mean_tit+FRE_se_tit), width=.8,col="blue") +
   geom_point(size = 2,col="blue")+
   
   labs(title= "Limnology",
        x = NULL, 
        y = "Readability [mean +/- S.E.]") +
   
   scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 
   
   theme_bw()+
   theme(
     legend.position = "none",
     axis.title = element_text(size = 10),
     axis.text.x = element_text(size = 10),
     axis.text.y = element_text(size = 10),
     panel.grid = element_blank(),
     plot.caption = element_text(size = 10, color = "gray50"),
     plot.title = element_text(face="bold", size=12)
   )
)

(FRE_t_OC <- ggplot(data=db_OC_year, aes(x=Publication.Year, y=FRE_mean_tit)) + 

   geom_line(aes(x=Publication.Year, y=FRE_mean_tit),linetype = 1,col="blue",alpha=1) + 
   
   geom_errorbar(aes(ymin=FRE_mean_tit-FRE_se_tit, ymax=FRE_mean_tit+FRE_se_tit), width=.8,col="blue") +
   geom_point(size = 2,col="blue")+
   
   labs(title= "Oceanography",
        x = NULL, 
        y = "Readability [mean +/- S.E.]") +
   
   scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 
   
   theme_bw()+
   theme(
     legend.position = "none",
     axis.title = element_text(size = 10),
     axis.text.x = element_text(size = 10),
     axis.text.y = element_text(size = 10),
     panel.grid = element_blank(),
     plot.caption = element_text(size = 10, color = "gray50"),
     plot.title = element_text(face="bold", size=12)
   )
)

(FRE_t_LA <- ggplot(data=db_LA_year, aes(x=Publication.Year, y=FRE_mean_tit)) + 

   geom_line(aes(x=Publication.Year, y=FRE_mean_tit),linetype = 1,col="blue",alpha=1) + 
   
   geom_errorbar(aes(ymin=FRE_mean_tit-FRE_se_tit, ymax=FRE_mean_tit+FRE_se_tit), width=.8,col="blue") +
   geom_point(size = 2,col="blue")+
   
   labs(title= "Lake Ecology",
        x = NULL, 
        y = "Readability [mean +/- S.E.]") +
   
   scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 
   
   theme_bw()+
   theme(
     legend.position = "none",
     axis.title = element_text(size = 10),
     axis.text.x = element_text(size = 10),
     axis.text.y = element_text(size = 10),
     panel.grid = element_blank(),
     plot.caption = element_text(size = 10, color = "gray50"),
     plot.title = element_text(face="bold", size=12)
   )
)

(FRE_t_HY <- ggplot(data=db_HY_year, aes(x=Publication.Year, y=FRE_mean_tit)) + 

   geom_line(aes(x=Publication.Year, y=FRE_mean_tit),linetype = 1,col="blue",alpha=1) + 
   
   geom_errorbar(aes(ymin=FRE_mean_tit-FRE_se_tit, ymax=FRE_mean_tit+FRE_se_tit), width=.8,col="blue") +
   geom_point(size = 2,col="blue")+
   
   labs(title= "Hydrobiology",
        x = NULL, 
        y = "Readability [mean +/- S.E.]") +
   
   scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 
   
   theme_bw()+
   theme(
     legend.position = "none",
     axis.title = element_text(size = 10),
     axis.text.x = element_text(size = 10),
     axis.text.y = element_text(size = 10),
     panel.grid = element_blank(),
     plot.caption = element_text(size = 10, color = "gray50"),
     plot.title = element_text(face="bold", size=12)
   )
)

# four plots together

FRE_t1 <- ggplotGrob(FRE_t_LI)
FRE_t2 <- ggplotGrob(FRE_t_OC)
FRE_t3 <- ggplotGrob(FRE_t_LA)
FRE_t4 <- ggplotGrob(FRE_t_HY)
grid.arrange(
	FRE_t1, FRE_t2, FRE_t3, FRE_t4,
	nrow=4,
	top="Readability of Title"
	)


# Number of authors --------------------------------------------------------

(Nau_LI <- ggplot(data=db_LI_year, aes(x=Publication.Year, y=N_aut_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_aut_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_aut_mean-N_aut_se, ymax=N_aut_mean+N_aut_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
        
    labs(title= "Limnology",
         x = NULL, 
         y = "Number of authors [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960, to=2020, by=10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 
    
    ylim(1,7) + 
                       
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

(Nau_OC <- ggplot(data=db_OC_year, aes(x=Publication.Year, y=N_aut_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_aut_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_aut_mean-N_aut_se, ymax=N_aut_mean+N_aut_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
        
    labs(title= "Oceanography",
         x = NULL, 
         y = "Number of authors [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 
    
    ylim(1,7) + 
    
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

(Nau_LA <- ggplot(data=db_LA_year, aes(x=Publication.Year, y=N_aut_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_aut_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_aut_mean-N_aut_se, ymax=N_aut_mean+N_aut_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
        
    labs(title= "Lake Ecology",
         x = NULL, 
         y = "Number of authors [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 
    
    ylim(1,7) + 
    
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

(Nau_HY <- ggplot(data=db_HY_year, aes(x=Publication.Year, y=N_aut_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_aut_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_aut_mean-N_aut_se, ymax=N_aut_mean+N_aut_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
        
    labs(title= "Hydrobiology",
         x = NULL, 
         y = "Number of authors [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 
    
    ylim(1,7) + 
    
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

# four plots together

Nau_1 <- ggplotGrob(Nau_LI)
Nau_2 <- ggplotGrob(Nau_OC)
Nau_3 <- ggplotGrob(Nau_LA)
Nau_4 <- ggplotGrob(Nau_HY)
grid.arrange(
	Nau_1, Nau_2, Nau_3, Nau_4,
	nrow=4,
	top="Number of authors"
	)


# Number of pages ---------------------------------------------------------

(Npa_LI <- ggplot(data=db_LI_year, aes(x=Publication.Year, y=N_pages_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_pages_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_pages_mean-N_pages_se, ymax=N_pages_mean+N_pages_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
        
    labs(title= "Limnology",
         x = NULL, 
         y = "Number of authors [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,30) + 
    
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

(Npa_OC <- ggplot(data=db_OC_year, aes(x=Publication.Year, y=N_pages_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_pages_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_pages_mean-N_pages_se, ymax=N_pages_mean+N_pages_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
            
    labs(title= "Oceanography",
         x = NULL, 
         y = "Number of authors [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,30) + 
    
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

(Npa_LA <- ggplot(data=db_LA_year, aes(x=Publication.Year, y=N_pages_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_pages_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_pages_mean-N_pages_se, ymax=N_pages_mean+N_pages_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
        
    labs(title= "Lake Ecology",
         x = NULL, 
         y = "Number of authors [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,30) + 
    
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

(Npa_HY <- ggplot(data=db_HY_year, aes(x=Publication.Year, y=N_pages_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_pages_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_pages_mean-N_pages_se, ymax=N_pages_mean+N_pages_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
        
    labs(title= "Hydrobiology",
         x = NULL, 
         y = "Number of authors [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,30) + 
    
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

# four plots together

Npa_1 <- ggplotGrob(Npa_LI)
Npa_2 <- ggplotGrob(Npa_OC)
Npa_3 <- ggplotGrob(Npa_LA)
Npa_4 <- ggplotGrob(Npa_HY)
grid.arrange(
	Npa_1, Npa_2, Npa_3, Npa_4,
	nrow=4,
	top="Number of pages"
	)


# Number of cited references ----------------------------------------------

(Ncr_LI <- ggplot(data=db_LI_year, aes(x=Publication.Year, y=N_citedref_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_citedref_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_citedref_mean-N_citedref_se, ymax=N_citedref_mean+N_citedref_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
    
    labs(title= "Limnology",
         x = NULL, 
         y = "Number of cited references [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,80) + 
    
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

(Ncr_OC <- ggplot(data=db_OC_year, aes(x=Publication.Year, y=N_citedref_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_citedref_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_citedref_mean-N_citedref_se, ymax=N_citedref_mean+N_citedref_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
    
    labs(title= "Oceanography",
         x = NULL, 
         y = "Number of cited references [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,80) + 
    
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

(Ncr_LA <- ggplot(data=db_LA_year, aes(x=Publication.Year, y=N_citedref_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_citedref_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_citedref_mean-N_citedref_se, ymax=N_citedref_mean+N_citedref_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
    
    labs(title= "Lake Ecology",
         x = NULL, 
         y = "Number of cited references [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,80) + 
    
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

(Ncr_HY <- ggplot(data=db_HY_year, aes(x=Publication.Year, y=N_citedref_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_citedref_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_citedref_mean-N_citedref_se, ymax=N_citedref_mean+N_citedref_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
    
    labs(title= "Hydrobiology",
         x = NULL, 
         y = "Number of cited references [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,80) + 
    
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

# four plots together

Ncr_1 <- ggplotGrob(Ncr_LI)
Ncr_2 <- ggplotGrob(Ncr_OC)
Ncr_3 <- ggplotGrob(Ncr_LA)
Ncr_4 <- ggplotGrob(Ncr_HY)
grid.arrange(
	Ncr_1, Ncr_2, Ncr_3, Ncr_4,
	nrow=4,
	top="Number of cited references"
	)


# Number of citations -----------------------------------------------------

(Nci_LI <- ggplot(data=db_LI_year, aes(x=Publication.Year, y=N_citations_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_citations_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_citations_mean-N_citations_se, ymax=N_citations_mean+N_citations_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
    
    labs(title= "Limnology",
         x = NULL, 
         y = "Number of citations [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,80) + 
    
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

(Nci_OC <- ggplot(data=db_OC_year, aes(x=Publication.Year, y=N_citations_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_citations_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_citations_mean-N_citations_se, ymax=N_citations_mean+N_citations_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
        
    labs(title= "Oceanography",
         x = NULL, 
         y = "Number of citations [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,80) + 
    
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

(Nci_LA <- ggplot(data=db_LA_year, aes(x=Publication.Year, y=N_citations_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_citations_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_citations_mean-N_citations_se, ymax=N_citations_mean+N_citations_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
    
    labs(title= "Lake Ecology",
         x = NULL, 
         y = "Number of citations [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,80) + 
    
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

(Nci_HY <- ggplot(data=db_HY_year, aes(x=Publication.Year, y=N_citations_mean)) + 
    
      geom_line(aes(x=Publication.Year, y=N_citations_mean),linetype = 1,col="blue",alpha=1) + 
      geom_errorbar(aes(ymin=N_citations_mean-N_citations_se, ymax=N_citations_mean+N_citations_se), width=.8,col="blue") +
    geom_point(size = 2,col="blue")+
    
    labs(title= "Hydrobiology",
         x = NULL, 
         y = "Number of citations [mean +/- S.E.]") + 
    
    scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,80) + 
    
    theme_bw()+
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 10, color = "gray50"),
      plot.title = element_text(face="bold", size=12)
    )
)

# four plots together

Nci_1 <- ggplotGrob(Nci_LI)
Nci_2 <- ggplotGrob(Nci_OC)
Nci_3 <- ggplotGrob(Nci_LA)
Nci_4 <- ggplotGrob(Nci_HY)
grid.arrange(
	Nci_1, Nci_2, Nci_3, Nci_4,
	nrow=4,
	top="Number of citations"
	)


# Number of papers --------------------------------------------------------

(Npa_LI <- ggplot(data=db_LI_year, aes(x=Publication.Year, y=N_papers)) + 

   geom_line(aes(x=Publication.Year, y= N_papers),linetype = 1,col="blue",alpha=1) + 
   
   #geom_errorbar(aes(ymin=FRE_mean_abs-FRE_se_abs, ymax=FRE_mean_abs+FRE_se_abs), width=.8,col="blue") +
   #geom_point(size = 2,col="blue")+
   
   labs(title= "Limnology",
        x = NULL, 
        y = "Number of papers") +
   
   scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,2000) + 
   
   theme_bw()+
   theme(
     legend.position = "none",
     axis.title = element_text(size = 10),
     axis.text.x = element_text(size = 10),
     axis.text.y = element_text(size = 10),
     panel.grid = element_blank(),
     plot.caption = element_text(size = 10, color = "gray50"),
     plot.title = element_text(face="bold", size=12)
   )
)

(Npa_OC <- ggplot(data=db_OC_year, aes(x=Publication.Year, y=N_papers)) + 

   geom_line(aes(x=Publication.Year, y= N_papers),linetype = 1,col="blue",alpha=1) + 
   
   #geom_errorbar(aes(ymin=FRE_mean_abs-FRE_se_abs, ymax=FRE_mean_abs+FRE_se_abs), width=.8,col="blue") +
   #geom_point(size = 2,col="blue")+
   
   labs(title= "Oceanography",
        x = NULL, 
        y = "Number of papers") +
   
   scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,2000) + 
   
   theme_bw()+
   theme(
     legend.position = "none",
     axis.title = element_text(size = 10),
     axis.text.x = element_text(size = 10),
     axis.text.y = element_text(size = 10),
     panel.grid = element_blank(),
     plot.caption = element_text(size = 10, color = "gray50"),
     plot.title = element_text(face="bold", size=12)
   )
)

(Npa_LA <- ggplot(data=db_LA_year, aes(x=Publication.Year, y=N_papers)) + 

   geom_line(aes(x=Publication.Year, y= N_papers),linetype = 1,col="blue",alpha=1) + 
   
   #geom_errorbar(aes(ymin=FRE_mean_abs-FRE_se_abs, ymax=FRE_mean_abs+FRE_se_abs), width=.8,col="blue") +
   #geom_point(size = 2,col="blue")+
   
   labs(title= "Lake Ecology",
        x = NULL, 
        y = "Number of papers") +
   
   scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,2000) + 
   
   theme_bw()+
   theme(
     legend.position = "none",
     axis.title = element_text(size = 10),
     axis.text.x = element_text(size = 10),
     axis.text.y = element_text(size = 10),
     panel.grid = element_blank(),
     plot.caption = element_text(size = 10, color = "gray50"),
     plot.title = element_text(face="bold", size=12)
   )
)

(Npa_HY <- ggplot(data=db_HY_year, aes(x=Publication.Year, y=N_papers)) + 

   geom_line(aes(x=Publication.Year, y= N_papers),linetype = 1,col="blue",alpha=1) + 
   
   #geom_errorbar(aes(ymin=FRE_mean_abs-FRE_se_abs, ymax=FRE_mean_abs+FRE_se_abs), width=.8,col="blue") +
   #geom_point(size = 2,col="blue")+
   
   labs(title= "Hydrobiology",
        x = NULL, 
        y = "Number of papers") +
   
   scale_x_continuous(breaks = c(seq(from=1960,to=2020,by = 10)), 
                       labels = c("1960","1970","1980","1990","2000","2010","2020"))+ 

    ylim(1,2000) + 
   
   theme_bw()+
   theme(
     legend.position = "none",
     axis.title = element_text(size = 10),
     axis.text.x = element_text(size = 10),
     axis.text.y = element_text(size = 10),
     panel.grid = element_blank(),
     plot.caption = element_text(size = 10, color = "gray50"),
     plot.title = element_text(face="bold", size=12)
   )
)

# four plots together

Npa_1 <- ggplotGrob(Npa_LI)
Npa_2 <- ggplotGrob(Npa_OC)
Npa_3 <- ggplotGrob(Npa_LA)
Npa_4 <- ggplotGrob(Npa_HY)
grid.arrange(
	Npa_1, Npa_2, Npa_3, Npa_4,
	nrow=4,
	top="Number of papers"
	)


# ANALYSES FOR NUMBER OF PAPERS

# prepare dataset
merged <- rbind(db_LI_year, db_OC_year, db_LA_year, db_HY_year)
merged$topic <- as.factor(c(rep('limnology',nrow(db_LI_year)),
				rep('oceanography',nrow(db_OC_year)),
				rep('lake ecology',nrow(db_LA_year)),
				rep('hydrobiology',nrow(db_HY_year))))
summary(merged)
names(merged)

# load libraries
library('car')
library('emmeans')
library('gam')
library('MASS')
library('performance')

# run gam model
m_pap <- gam::gam(N_papers ~ topic*s(Publication.Year),
			 family=poisson, data=merged)
performance::check_model(m_pap)
summary(m_pap)
pairs(emmeans::emmeans(m_pap, ~topic*s(Publication.Year)),
	simple="topic")

# run nb.glm model	
m_pap2 <- MASS::glm.nb(N_papers ~ topic*Publication.Year, data=merged)
performance::check_model(m_pap2)
car::Anova(m_pap2)
pairs(emmeans::emmeans(m_pap2, ~topic*Publication.Year),
	simple="topic")

# plot the results
ggplot2::ggplot(merged, aes(x=Publication.Year, y=N_papers)) + 
	theme_light() +
	geom_point(aes(colour=topic)) + 
	geom_smooth(aes(colour=topic)) +
	scale_y_continuous(trans='log2') +
	scale_x_continuous(breaks = c(seq(from=1960, to=2020, by=10)), 
         labels=c("1960","1970","1980","1990","2000","2010","2020")) +
    labs(x=NULL, y="Number of papers") +
    scale_color_brewer(palette="Paired")