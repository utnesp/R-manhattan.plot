# R-manhattan.plot
Create a nice looking manhattan plot of significantly up- and downregulated genes. 

## INSTALLATION

```R
devtools::install_github("utnesp/Easy-bioMart")
devtools::install_github("utnesp/R-manhattan.plot")
library(easybiomart)
library(plot.manhattan)
```

## USAGE

```R
> counts[1:6,1:6]
                sample_1 sample_2 sample_3 sample_4 sample_5 sample_6
ENSG00000227232      109      117       63      179       67       81
ENSG00000268903       18        3        4       16        1       18
ENSG00000269981       14        1        5       12        1        7
ENSG00000279457      463      472      261      247      350      475
ENSG00000228463        0       14       12       65       21        8
ENSG00000237094      135       47       30       32       47       99

> head(clinical.attributes)
         Eventfree Alive           MNA
sample_1     EVENT  DEAD NON-AMPLIFIED
sample_2     EVENT  DEAD     AMPLIFIED
sample_3     EVENT  DEAD     AMPLIFIED
sample_4     EVENT  DEAD     AMPLIFIED
sample_5     EVENT  DEAD     AMPLIFIED
sample_6     EVENT  DEAD     AMPLIFIED
```

If we would like to create a manhattan plot based whether or not the tumor sample is MYCN amplified or not, we would set which.attribute = 3. We also need to specifiy where to save the output with file.loc. The fold change, P and FDR values will be based on the output from edgeR.

```R
plot.manhattan(counts, clinical.attributes, file.loc = "/Volumes/localdisc/output.folder/")
```

We can also specify other options such as:

```R
open.pdf = T,open.html = T          ## will open the manhattan plot after the function has completed 
reduce.point.size.by = 2            ## useful when point sizes are to0 big or too low
highligt.top.genes = T, nr_top_genes = 20, nr_bottom_genes = 20           ## for labelling xx genes of interrest
which.biotypes = "protein_coding"   ## if interested in comparing only protein_coding genes
which.biotypes = "protein_coding", inverse.biotypes = F                   ## if interested in only long non-coding RNAs
reverseGroup = T                    ## if the comparison made by edgeR is not what we wanted, e.g. we would like AMPLIFIED / NON-AMPLIFIED
instead of NON-AMPLIFIED / AMPLIFIED
name = ""                           ## if you would like to specify the output names of each file manually. Defaults to column name of clinical attribute
title.name = ""                     ## Title name in plot. Defaults to column name of clinical attribute
count.mean.cutoff = 10              ## Cutoff before running edgeR analysis
logFC.cutoff = 0.8                  ## Cutoff before correcting p values
Pvalue.cutoff = 0.05                ## Plot only genes that satisfy this criteria
write.all.genes = T                 ## Wheter to write all genes in edgeR output, or to only write those below Pvalue.cutoff
dpi = 300, width = 297, height = 210, units = "mm"    ##  pixels, width, height of pdf
order.by.pvalue = F                 ## order and plot by pvalue rather than FDR
write.cpm = F                       ## write count matrix to output folder

```


The files created by this function will be saved in the folder specified by you in file.loc.

Example showing up- and downregulated genes with FDR < 0.05 in MYCN amplified vs MYCN non-amplified samples. 
The point size correspons to how significant the genes are. In this case we see that MYCN is signifcantly upregulated in the comparison AMPLIFIED vs NON-AMPLFIED tumors. The different point types will also correspond to whether or not the gene of interest is a protein coding gene or not (circle = protein_coding, triangle = non-coding)


![alt text](https://github.com/utnesp/R-manhattan.plot/blob/master/MNA.labeled.jpg "Example showing up- and downregulated genes with FDR < 0.05")



The function also creates an interactive html plot (using ggplotly) as well as the output from edgeR. 

Depends: 
edgeR, ggplot and ggplotly, and
https://github.com/utnesp/Easy-bioMart
