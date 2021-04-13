# Lake Varese files description

The Lake Varese is part of the AQVST project and consists of a metagenomic dataset sampled from organic matter, microplastic and water from Lago Varese and nearbies.
The annotation of the whole dataset can be found on [MetaStorm](http://bench.cs.vt.edu/MetaStorm/) by browsing the "Varese Lake microplastics" project.

- The files "LakeVarese_MetaPhlAn_family_scaling.csv" and "LakeVarese_MetaPhlAn_genus_scaling.csv" are the pivot table coming from the taxonomy annoation of the dataset against the MetaPhlAn database (v2) on MetaStorm.
- The file "pivot_filtering.Rmd" is the R script which has been used to filter Virus and Archaea taxa out of the taxonomy pivot tables. 
- The files "LakeVarese_MetaPhlAn_filtered_family_scaling.xlsx" and "LakeVarese_MetaPhlAn_filtered_genus_scaling.xlsx" are the output of "pivot_filtering.Rmd".
