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

Make the figure


```r
##pi
dist = 5000
abc = c('C','C','A','B')

par(mfrow=c(3,1), mar = c(5,5,3,1))
for (i in c(3,4,2)) {
mysnp = myhits[i,]
#print(mysnp)
mydiv = dplyr::filter(div, scaf == mysnp$scaf)
plot(mydiv$midpoint, mydiv$pi,xlim = c(mysnp$pos-dist,mysnp$pos+dist), ylab = "", col = "darkgray", xaxt="n", cex=1.5, yaxt="n", ylim=c(0,0.06), bty="n", xlab = 'kB', lwd=2)
#abline(v=mysnp$ps, col = mysnp$dispColor, lwd=2)
abline(v=mysnp$pos, col = 'black', lwd=2, lty=2)
if(i == 2){abline(v = myhits[1,]$pos, col = "black", lwd=2, lty=2)}
abline(h=pi.975, col = "darkgrey", lwd=2, lty=2)
axis(1, at=seq(mysnp$pos-dist,mysnp$pos+dist,250), labels = seq(-dist,dist,by=250)/1000, cex.axis=1.5, padj=1)
#axis(1, at = seq(mysnp$ps-dist,mysnp$ps+dist,1000), labels = seq(mysnp$ps-dist,mysnp$ps+dist,1000)/1000)
axis(2, las=2, cex.axis=1.5, at = c(0:2)/20)
mtext(expression(pi),side=2, cex=1, line=4)
mtext(abc[i], side=3, adj=-0.07, cex=1.3, line=1)

##add in the second eqtl
}
```

![](popgen-of-eQTLs_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
#tajd
##pi
dist = 5000
abc = c('C','C','A','B')

par(mfrow=c(3,1), mar = c(5,7,3,1))
for (i in c(3,4,2)) {
mysnp = myhits[i,]
#print(mysnp)
mydiv = dplyr::filter(div, scaf == mysnp$scaf)
plot(mydiv$midpoint, mydiv$tajd,xlim = c(mysnp$pos-dist,mysnp$pos+dist), ylab = "", col = "darkgray", xaxt="n", cex=1.5, yaxt="n", ylim=c(-2.6,2), bty="n", xlab = 'kB', lwd=2)
#abline(v=mysnp$ps, col = mysnp$dispColor, lwd=2)
abline(v=mysnp$pos, col = 'black', lwd=2, lty=2)
if(i == 2){abline(v = myhits[1,]$pos, col = "black", lwd=2, lty=2)}
abline(h=pi.975, col = "darkgrey", lwd=2, lty=2)
axis(1, at=seq(mysnp$pos-dist,mysnp$pos+dist,250), labels = seq(-dist,dist,by=250)/1000, cex.axis=1.5, padj=1)
#axis(1, at = seq(mysnp$ps-dist,mysnp$ps+dist,1000), labels = seq(mysnp$ps-dist,mysnp$ps+dist,1000)/1000)
axis(2, las=2, cex.axis=1.5, at = c(-2:2))
mtext("Tajima's D",side=2, cex=1, line=4)
mtext(abc[i], side=3, adj=-0.07, cex=1.3, line=1)

##add in the second eqtl
}
```

![](popgen-of-eQTLs_files/figure-html/unnamed-chunk-2-2.png)<!-- -->

