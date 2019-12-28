---
title: "matrix-eQTL"
author: "em"
date: "December 28, 2019"
output:
  html_document:
    keep_md: yes
---



##Get the batch-corrected expression data together for mapping

```r
load('data/wgcna-files/combatOutput.rda')

mycbnames = read.table("data/wgcna-files/mycb-names", stringsAsFactors = F)

qn <- function(x){qqnorm(x, plot.it=FALSE)$x}
# nnormed = data.frame(t(apply(outtable[,-1], 1,qn)))
normcb = apply(mycb, 1, qn)#quantile normalize, now columns genes, rows are inds


mycbdf = data.frame(id = mycbnames$x, normcb, stringsAsFactors = F)

matrixeqtlids = data.frame(V1 = strsplit('10    100    101    102_long    105_clipinplace_sorted    106s_clipped    107    108_3    109    110    111    112    113    114_long    115    117s_clipped    118    11s_clipped    121    123    124_3    125_clipped.bam    126s_clipped    128    129s_clipped    13    131s_clipped    132    133    135    136    137s_clipped    138    139    140    141    142    143s_clipped    144s_clipped    147    148    14_clipped.bam    15    151    152    153    154    155_clipped.bam    156    157    158    160    161_clipped.bam    162    163_clipped    165_clipped.bam    167    16_clipped    17    170_clipped    174_clipped.bam    175_clipped    176    177    178    18    181    182_3    183    184_3    186    187    189    190    192    193    194    195    197    198    1_clipped    20    200    202_long    203_long    204_long    23s_clipped    24    25    26    27    28    29_clipped    30    31    33    34    35    36    38    39s_clipped    41    43s_clipped    44s_clipped    47    49_clipped.bam    4s_clipped    50_long    52    54_clipped    55_clipped    58    5_clipped    61    63    64_replacement    65_clipped    66_clipped    67    70s_clipped    71    72    74    75    76    78s_clipped    79    7s_clipped    8    80s_clipped    81s_clipped    82s_clipped    83    85_long    86    89    9    91    92    93    94    95    96    97_3    98_clipped    99_clipped', split = "    ")[[1]], stringsAsFactors = F)



meExp = data.frame(c('id', row.names(mycb)),t(dplyr::left_join(matrixeqtlids, mycbdf, by = c('V1'='id'))), stringsAsFactors = F)
write.table(meExp, file = 'data/matrixeqtl-files/cbExp', row.names=F, quote=F, col.names=F)
```


