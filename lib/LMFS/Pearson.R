getwd()
library(ggpubr)

####single sample (2 functions, same results)
tot2<-read.csv2("ARG-MRG2.csv")
tot2
test<-cor(tot2$ARGs, tot2$MRGs, method = "pearson")
test
test2<-cor.test(tot2$ARGs, tot2$MRGs, method = "pearson")
test2

ggscatter(tot2, x = "ARGs", y = "MRGs",
          add = "reg.line", cor.coef = TRUE, cor.method = "pearson",
          xlab = "ARG relative total abundance (normalized to 16S rRNA)", ylab = "MRG relative total abundance (normalized to 16S rRNA)")

####averaged samples
tot3<-aggregate(tot2[,c(2,3)], list(tot2$Sample), mean)
tot3
test3<-cor(tot3$ARGs, tot3$MRGs, method = "pearson")
test3
test4<-cor.test(tot3$ARGs, tot3$MRGs, method = "pearson")
test4

ggscatter(tot3, x = "ARGs", y = "MRGs",
          add = "reg.line", cor.coef = TRUE, cor.method = "pearson",
          xlab = "ARG relative total abundance (normalized to 16S rRNA)", ylab = "MRG relative total abundance (normalized to 16S rRNA)")
