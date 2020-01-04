---
title: "popgen"
author: "em"
date: "January 4, 2020"
output:
  html_document:
    keep_md: yes
---




```r
##read in the diversity data
div = read.table("data/scall_014diversity_final", stringsAsFactors = F, header=T)
names(div) = c("scaf","midpoint",'pi',	'theta',	'tajd',	'divergence')

## get upper tail of 95% CI
pi.975 = quantile(div$pi, probs = c(0.975))
d.975 = quantile(div$tajd, probs = c(0.975), na.rm=T)


## read in the eqlts
ceqtls = read.table('Table1', header=T, stringsAsFactors = F)
ceqtls$scaf = sapply(ceqtls$loc, function(x){ strsplit(as.character(x), split=":")[[1]][1]  })
ceqtls$pos = sapply(ceqtls$loc, function(x){ as.numeric(strsplit(as.character(x), split=":")[[1]][2])  })
ceqtls$myIndex = as.numeric(row.names(ceqtls))

#make genomic ranges object for coexpression eQTLs
eqtlranges = GRanges(seqname = ceqtls$scaf, ranges = IRanges(start=ceqtls$pos, width=1))

#make genomic ranges for div data
myranges = GRanges(seqname=div$scaf, ranges = IRanges(start=div$midpoint-250, end=div$midpoint+249))

##find overlaps
myOverlaps = as.matrix(findOverlaps(myranges,eqtlranges))
div$myIndex = as.numeric(row.names(div))
myhits = dplyr::inner_join(as.data.frame(myOverlaps), ceqtls, by = c("subjectHits" = "myIndex")) 
myhitsdiv = dplyr::inner_join(myhits, div, by=c("queryHits"="myIndex"))

## information about pi and tajd around coexpression eQTLs
summary(myhitsdiv$pi)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.006620 0.006707 0.008884 0.011612 0.013789 0.022058
```

```r
summary(myhitsdiv$d)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0000  0.0000  0.1183  0.1809  0.2992  0.4872
```
