# R-manhattan.plot
Create a manhattan plot using ggplot and ggplotly. 

```R
> counts[1:5,1:5]
                sample_1 sample_2 sample_3 sample_4 sample_5
ENSG00000227232      109      117       63      179       67
ENSG00000268903       18        3        4       16        1
ENSG00000269981       14        1        5       12        1
ENSG00000279457      463      472      261      247      350
ENSG00000228463        0       14       12       65       21

> head(clinical.attributes)
         Eventfree Alive           MNA
sample_1     EVENT  DEAD NON-AMPLIFIED
sample_2     EVENT  DEAD     AMPLIFIED
sample_3     EVENT  DEAD     AMPLIFIED
sample_4     EVENT  DEAD     AMPLIFIED
sample_5     EVENT  DEAD     AMPLIFIED
sample_6     EVENT  DEAD     AMPLIFIED

# which.attribute correspons in this case to column named MNA
# if the points are to big, we can reduce the size using reduce.point.size.by
# we can use highlight.top.genes to label the most significantly up- and downregulated genes
# we can specify whether or not to highlight protein_coding genes, or all genes except protein_coding genes
# the differential expression depends on edgeR. If the comparison set should be reversed, use reverseGroup. Ex. if you would like AMPLIFIED / NON-AMPLIFIED rather than NON-AMPLIFIED / AMPLIFIED
plot.manhattan(counts, clinical.attributes, file.loc = "/Volumes/localdisc/output.folder/", open.pdf = F, which.attribute = 3, reduce.point.size.by = 2, highligt.top.genes = T, nr_top_genes = 20, nr_bottom_genes = 20, which.biotypes = "protein_coding", inverse.biotypes = F, reverseGroup = T)

```
The files created by this function will be saved in the folder specified by you in file.loc.

Example showing up- and downregulated genes with FDR < 0.05 in MYCN amplified vs MYCN non-amplified samples. 
The point size correspons to how significant the genes are. In this case we see that MYCN is signifcantly upregulated in the comparison AMPLIFIED vs NON-AMPLFIED tumors
![alt text](https://github.com/utnesp/R-manhattan.plot/blob/master/MNA.labeled.pdf "Example showing up- and downregulated genes with FDR < 0.05")

The function also creates an interactive html plot (using ggplotly) and the output from edgeR. 
