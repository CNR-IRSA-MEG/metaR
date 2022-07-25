setwd("~/Dropbox/R scripts/Black_sea/int_genes")

library(readr)
library(dplyr)
library(stringr)



## read the outputs,  repeat or loop the code
SAMPLE_ID_int <- read_delim("SAMPLE_ID_int.txt", delim = "|", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE) 
colnames(SAMPLE_ID_int) <- c('gnl','code1','database','code2','isolate','gene') # rename cols
SAMPLE_ID_int_filt <- SAMPLE_ID_int %>% filter(grepl('int|Int', gene))  # select rows containing int or Int
SAMPLE_ID_int_filt_summary <- SAMPLE_ID_int_filt %>% group_by(gene) %>% summarise(count=n()) # count int-like groups hits
colnames(SAMPLE_ID_int_filt_summary) <- c('gene','SAMPLE_ID') # rename columns as per samples


## now i merge tables

int_merged  = SAMPLE_ID_int_filt_summary %>% 
  full_join(SAMPLE_ID2_int_filt_summary, by = 'gene') %>%
  full_join(SAMPLE_ID3_int_filt_summary, by = 'gene') %>%
  full_join(SAMPLE_ID4_int_filt_summary, by = 'gene')  # and so on
 

# i replace NAs with zeros
int_merged[is.na(int_merged)] <- 0 # NAs as 0
View(int_merged)
# 
# save file
write.table(int_merged, file = "int_merged.csv", sep = ",", col.names = NA, qmethod = "double")


#####
