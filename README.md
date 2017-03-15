# R-manhattan.plot
Create a manhattan plot using edgeR, ggplot and ggplotly. 

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
open.pdf = F                        ## will open the manhattan plot after the function has completed 
reduce.point.size.by = 2            ## useful when point sizes are to0 big or too low
highligt.top.genes = T, nr_top_genes = 20, nr_bottom_genes = 20           ## for labelling xx genes of interrest
which.biotypes = "protein_coding"   ## if interested in comparing only protein_coding genes
which.biotypes = "protein_coding", inverse.biotypes = F                   ## if interested in only long non-coding RNAs
reverseGroup = T                    ## if the comparison made by edgeR is not what we wanted, e.g. we would like AMPLIFIED / NON-AMPLIFIED
instead of NON-AMPLIFIED / AMPLIFIED
```


The files created by this function will be saved in the folder specified by you in file.loc.

Example showing up- and downregulated genes with FDR < 0.05 in MYCN amplified vs MYCN non-amplified samples. 
The point size correspons to how significant the genes are. In this case we see that MYCN is signifcantly upregulated in the comparison AMPLIFIED vs NON-AMPLFIED tumors
![alt text](https://github.com/utnesp/R-manhattan.plot/blob/master/MNA.labeled.pdf "Example showing up- and downregulated genes with FDR < 0.05")

The function also creates an interactive html plot (using ggplotly) and the output from edgeR. 
