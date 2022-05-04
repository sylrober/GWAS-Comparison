# GWAS-Comparison
R code to compare GWAS data of Alzheimer's and Depression

Description: First, the script imports the depression data set, which I titled dep.txt.
The data gets sorted by p value and thresholds are established. Rsid is isolated for each threshold. Then, it imports and does the same with the Alzheimer’s data set. Since rsid is missing, it creates csv files of each threshold that were used to query for rs id on Kaviar- http://db.systemsbiology.net/kaviar/cgi-pub/Kaviar.pl. The results were saved as text files and imported for the venn digram of rsid at the 3 thresholds. 

Next, genes were obtained from the snp’s. Overlapping genes were then obtained. A vector was created of the intersection.

All the significant GO and reactome results from enrichr were then used as keys in bioconductor to get gene symbols, and the gene overlap was obtained to connect genes to pathways. Then, incidence matrices were made from these results in excel. 

The incidence matrices were imported as csv files for heatmaps visualization, and then igraph objects were created and connected to Cytoscape.
A Fisher’s exact test was then added to find the statistical significance of the 13 overlapping genes analyzed. 

Required files: 
Depression dataset:
https://datashare.ed.ac.uk/bitstream/handle/10283/3203/PGC_UKB_depression_genome-wide.txt?sequence=3&isAllowed=y

Alzheimer’s dataset:
https://ctg.cncr.nl/documents/p1651/PGCALZ2sumstatsExcluding23andMe.txt.gz

-alz.txt, alz2.txt, alz3.txt – rs id for Alzheimer’s thresholds 
-gocell.csv, gogenes1.csv, gomol.csv, re.csv – incidence matrices made to create heatmaps and networks

Required packages:
tidyverse, dplyr - to filter

VennDiagram, gplots - to compare overlap

BiocManager, biomaRt - for identifying genes

GO.db, org.Hs.eg.db - for genes from GO pathways

ReactomeContentService4R- for genes from Reactome pathways

igraph - for creating networks

RCy3 - to transfer networks to Cytoscape
*Cytoscape software is required to be opened for RCy3 to work

GeneOverlap - used for Fisher’s exact test

Execution:
1)    Open R studio
2)    Download required files. 
3)    Change working directory to location of files
4)    Open script file and run script
*Since the data files are large, the script takes a while to run. It’s best to run it in small chunks at a time. 

Output files:
Thresholds of Alzheimer’s data used to obtain rs id from Kaviar:
-t1.csv
-t2.csv
-t3.csv

6 venn diagrams output in R (3 of SNP overlap and 3 of gene overlap)

4 heat maps:  
-gobio.png
-gomol.png
-gocell.png
-re.png

4 Cytoscape networks:
-renetwork (reactome)
-go1network (biological)
-go2network (molecular)
-go3network (cellular)
