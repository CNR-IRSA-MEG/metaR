# Black Sea files description

The Black Sea dataset is a collection of aquatic metagenomic read files from the NCBI projects PRJNA638805 and PRJNA649215.
The annotation of the whole dataset can be found on [MetaStorm](http://bench.cs.vt.edu/MetaStorm/) by browsing the "Black Sea resistome" project.

- The file "pivot_creation.R" is an R script which was used to merge all annotated MAGs (DIAMOND BLASTP `.m8` tabular output), both against the BacMet and the DeepARG database, into two separate pivot table, one for the antimicrobial resistance functions and the other for the metal resistance functions.
Two additional pivots per kind were created, one reporting counts of the annotated functions and the other reporting normalized counts as scaling (values between 0 and 100).
