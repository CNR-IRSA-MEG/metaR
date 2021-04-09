#!/usr/bin/Rscript

#Loading libraries quietly
#If you need functions from both plyr and dplyr, please load plyr first, then dplyr
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(purrr)) # For map()
suppressPackageStartupMessages(library(stringr)) # For str_remove()
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(openxlsx)) # For writing to .xlsx file

#setting current directory as constant
BASE_GLOB = '.'

# Loads default BLOSUM62 alignment table (from BLAST, LAST, DIAMOND, etc.)
# Reads all .argdb.m8 blast/diamond output into a long table
long_arg <- tibble(fname = Sys.glob("*.argdb.m8")) %>%
  mutate(
    d = map(
      fname,
      function(f) read_tsv(f, col_names = c(
        "contig_ID", "db_descriptor", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", 
        "evalue", "bitscore", "qcovhsp"
      ), col_types = cols(.default = col_double(), contig_ID = col_character(), db_descriptor = col_character())
      ) %>%
        mutate(sample = str_remove(f, '\\.cdhit\\.argdb\\.m8'))
    )
  ) %>%
  unnest(d) %>%
  select(-fname) %>%
  separate(db_descriptor, c("db_ID", "features", "db", "phenotype", "ARG"), sep = "\\|", fill = 'right', extra = 'drop') %>%
  filter(qcovhsp >= 70)

# Do the same for .bacmet.m8 files
long_mrg <- tibble(fname = Sys.glob("*.bacmet.m8")) %>%
  mutate(
    d = map(
      fname,
      function(f) read_tsv(f, col_names = c(
        "contig_ID", "db_descriptor", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", 
        "evalue", "bitscore", "qcovhsp"
      ), col_types = cols(.default = col_double(), contig_ID = col_character(), db_descriptor = col_character())
      ) %>%
        mutate(sample = str_remove(f, '\\.cdhit\\.bacmet\\.m8'))
    )
  ) %>%
  unnest(d) %>%
  select(-fname) %>%
  separate(db_descriptor, c("db_ID", "MRG", "features", "secondary_id", "phenotype"), sep = "\\|", fill = 'right', extra = 'drop') %>%
  filter(qcovhsp >= 70)

antibiotics_counts <- long_arg %>%
  select(sample, ARG)%>%  
  # Count the number of rows for combinations of sample and ARG. This introduces a new column: n.
  count(sample, ARG) %>%
  # Make it wide: genes as rows, samples as columns
  pivot_wider(
    names_from = sample,
    values_from = n,
    values_fill = 0
  ) %>%
  write.xlsx('black_sea_ARGs_counts.xlsx')

antibiotics_scaling <- long_arg %>%
  count(sample, ARG) %>%
  group_by(sample) %>%
  mutate(relab = n/sum(n) * 100) %>%
  ungroup() %>%
# Make it wide: genes as rows, samples as columns
  pivot_wider( -n,
    names_from = sample,
    values_from = relab,
    values_fill = 0
  ) %>%
  write.xlsx('black_sea_ARGs_scaling.xlsx')

metals_counts <- long_mrg %>%
  select(sample, MRG)%>%  
  # Count the number of rows for combinations of sample and ARG. This introduces a new column: n.
  count(sample, MRG) %>%
  # Make it wide: genes as rows, samples as columns
  pivot_wider(
    names_from = sample,
    values_from = n,
    values_fill = 0
  ) %>%
  write.xlsx('black_sea_MRGs_counts.xlsx')

metals_scaling <- long_mrg %>%
  count(sample, MRG) %>%
  group_by(sample) %>%
  mutate(relab = n/sum(n) * 100) %>%
  ungroup() %>%
  # Make it wide: genes as rows, samples as columns
  pivot_wider( -n,
               names_from = sample,
               values_from = relab,
               values_fill = 0
  ) %>%
  write.xlsx('black_sea_MRGs_scaling.xlsx')