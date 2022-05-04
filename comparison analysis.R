#import data
mydata=read.table('dep.txt', header = TRUE, sep = " ")

##Depression Genes
sorteddata=mydata[order(mydata$P),]
sorteddata

library(tidyverse)
library(dplyr)

#number of top SNPS 
threshold1 = sorteddata %>% filter(P <= 0.00000005)
threshold2=sorteddata %>% filter(P <= 0.0000005)
threshold3=sorteddata %>% filter(P <= 0.000005)

#rsid
dep1=threshold1$MarkerName
dep2=threshold2$MarkerName
dep3=threshold3$MarkerName

##Alzheimer's Genes
#number of top SNPs of alzheimer's data
mydata2=read.table('PGCALZ2sumstatsExcluding23andMe.txt', header = TRUE, sep = "\t")
sorteddata2=mydata2[order(mydata2$p),]
sorteddata2

#number of top SNPS 
t1 = sorteddata2 %>% filter(p <= 0.00000005)
t2=sorteddata2 %>% filter(p <= 0.0000005)
t3=sorteddata2 %>% filter(p <= 0.000005)

#save csv to use for website that finds rsid
write.csv(t1,'t1.csv')
write.csv(t2,'t2.csv')
write.csv(t3,'t3.csv')

#import results
alz=read.table('alz.txt',sep="\n")
as.data.frame(alz)
alz=alz[[1]]

alz2=read.table('alz2.txt',sep="\n")
as.data.frame(alz2)
alz2=alz2[[1]]

alz3=read.table('alz3.txt',sep="\n")
as.data.frame(alz3)
alz3=alz3[[1]]

library(VennDiagram)
#diagrams with color showing overlapping SNP's
venn.diagram(x = list(alz, dep1) ,
             category.names = c("GWAS Alzheimer's", "GWAS Depression"),
             filename = 't1.png',
             output=TRUE,
             imagetype="png", 
             scaled = FALSE,
             col = "black",
             fill = c("orange","blue"),
             cat.col = c("orange","blue"),
             cat.cex = 1,
             margin = 0.15,
             cat.pos=c(0,10),
             main= expression(Threshold ~ 5*e^-8))

venn.diagram(x = list(alz2, dep2) ,
             category.names = c("GWAS Alzheimer's", "GWAS Depression"),
             filename = 't2.png',
             output=TRUE,
             imagetype="png", 
             scaled = FALSE,
             col = "black",
             fill = c("orange","blue"),
             cat.col = c("orange","blue"),
             cat.cex = 1,
             margin = 0.15,
             cat.pos=c(0,10),
             main= expression(Threshold ~ 5*e^-7))

venn.diagram(x = list(alz3, dep3) ,
             category.names = c("GWAS Alzheimer's", "GWAS Depression"),
             filename = 't3.png',
             output=TRUE,
             imagetype="png", 
             scaled = FALSE,
             col = "black",
             fill = c("orange","blue"),
             cat.col = c("orange","blue"),
             cat.cex = 1,
             margin = 0.15,
             cat.pos=c(0,10),
             main= expression(Threshold ~ 5*e^-6))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(biomaRt)
#get genes from rsid
ensembl <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
mart <- useEnsembl(biomart = "ENSEMBL_MART_SNP", 
                   dataset = "hsapiens_snp")
depgenes1=getBM(attributes=c( "ensembl_gene_stable_id"), filters="snp_filter", values=dep1, mart=ensembl, uniqueRows=TRUE)
depgenes2=getBM(attributes=c( "ensembl_gene_stable_id"), filters="snp_filter", values=dep2, mart=ensembl, uniqueRows=TRUE)
depgenes3=getBM(attributes=c( "ensembl_gene_stable_id"), filters="snp_filter", values=dep3, mart=ensembl, uniqueRows=TRUE)

Agenes1=getBM(attributes=c( "ensembl_gene_stable_id"), filters="snp_filter", values=alz, mart=ensembl, uniqueRows=TRUE)
Agenes2=getBM(attributes=c( "ensembl_gene_stable_id"), filters="snp_filter", values=alz2, mart=ensembl, uniqueRows=TRUE)

Agenes3=getBM(attributes=c( "ensembl_gene_stable_id"), filters="snp_filter", values=alz3, mart=ensembl, uniqueRows=TRUE)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
m=getBM(attributes = c("hgnc_symbol"), values=depgenes1,filters = "ensembl_gene_id",mart = ensembl,uniqueRows=TRUE)
m2=getBM(attributes = c("hgnc_symbol"), values=depgenes2,filters = "ensembl_gene_id",mart = ensembl,uniqueRows=TRUE)
m3=getBM(attributes = c("hgnc_symbol"), values=depbiogenes3,filters = "ensembl_gene_id",mart = ensembl,uniqueRows=TRUE)
a=getBM(attributes = c("hgnc_symbol"), values=Agenes1,filters = "ensembl_gene_id",mart = ensembl,uniqueRows=TRUE)
a2=getBM(attributes = c("hgnc_symbol"), values=Agenes2,filters = "ensembl_gene_id",mart = ensembl,uniqueRows=TRUE)
a3=getBM(attributes = c("hgnc_symbol"), values=Agenes3,filters = "ensembl_gene_id",mart = ensembl,uniqueRows=TRUE)

#compare genes and find overlap for enrichr
library(gplots)

table1=venn(list(m[[1]],a[[1]]))
print(table1)

table2=venn(list(m2[[1]],a2[[1]]))
print(table2)

table3=venn(list(m3[[1]],a3[[1]]))
print(table3)

#save overlapping genes of this threshold as vector
genes=attr(table3,"intersections")$`A:B`

#get genes for go results
library(GO.db)
library(org.Hs.eg.db)
#Go BIO
results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0060333'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols1 <- unique(results$SYMBOL)
print(gene_symbols1)
results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0019886'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols2 <- unique(results$SYMBOL)
print(gene_symbols2)
results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0002495'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols3 <- unique(results$SYMBOL)
print(gene_symbols3)
results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0002478'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols4 <- unique(results$SYMBOL)
print(gene_symbols4)
results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0035610'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols5 <- unique(results$SYMBOL)
print(gene_symbols5)
results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0071346'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols6 <- unique(results$SYMBOL)
print(gene_symbols6)
results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0016926'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols7 <- unique(results$SYMBOL)
print(gene_symbols7)
results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0050852'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols8 <- unique(results$SYMBOL)
print(gene_symbols8)
results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0035608'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols9 <- unique(results$SYMBOL)
print(gene_symbols9)
results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0050851'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols10 <- unique(results$SYMBOL)
print(gene_symbols10)
results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0070585'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols11 <- unique(results$SYMBOL)
print(gene_symbols11)

a=venn(list(e=	genes,	b=	gene_symbols1	))
print(a)
b=venn(list(e=	genes,	b=	gene_symbols2	))
print(b)
c=venn(list(e=	genes,	b=	gene_symbols3	))
print(c)
d=venn(list(e=	genes,	b=	gene_symbols4	))
print(d)
e=venn(list(e=	genes,	b=	gene_symbols5	))
print(e)
f=venn(list(e=  genes,	b=	gene_symbols6	))
print(f)
g=venn(list(e=	genes,	b=	gene_symbols7	))
print(g)
h=venn(list(e=	genes,	b=	gene_symbols8	))
print(h)
i=venn(list(e=	genes,	b=	gene_symbols9	))
print(i)
j=venn(list(e=	genes,	b=	gene_symbols10	))
print(j)
k=venn(list(e=	genes,	b=	gene_symbols11	))
print(k)

gogenes1=read.csv('gogenes1.csv',header = TRUE)
row.names(gogenes1)=c('antigen processing & presentation of exogenous peptide antigen via MHC class II',
'antigen processing & presentation of peptide antigen via MHC class II',
'antigen processing & presentation of exogenous peptide antigen',
'protein side chain deglutamylation',
'protein desumoylation',
'T cell receptor signaling pathway',
'protein deglutamylation',
'antigen receptor-mediated signaling pathway',
'protein localization to mitochondrion') 

gogenes1[is.na(gogenes1)] <- 0

#Go MOLECULAR
results <- AnnotationDbi::select(org.Hs.eg.db, keys=c(' GO:0032395'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols1 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0016929'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols2 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0033130'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols3 <- unique(results$SYMBOL)

a=venn(list(e=	genes,	b=	gene_symbols1	))
print(a)
b=venn(list(e=	genes,	b=	gene_symbols2	))
print(b)
c=venn(list(e=	genes,	b=	gene_symbols3	))
print(c)

gogenes2=read.csv('gomol.csv',header = TRUE)
row.names(gogenes2)=c('SUMO-specific protease activity',
'acetylcholine receptor binding') 

#Go CELLULAR
results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0042613'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols1 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0042611'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols2 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0098553'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols3 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0071556'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols4 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0012507'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols5 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0030662'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols6 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0030658'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols7 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0098852'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols8 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0030669'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols9 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0030134'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols10 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0005765'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols11 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0045334'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols12 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0030665'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols13 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0032588'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols14 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0005764'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols15 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0030176'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols16 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0030666'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols17 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0030139'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols18 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0005802'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols19 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0030659'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols20 <- unique(results$SYMBOL)

results <- AnnotationDbi::select(org.Hs.eg.db, keys=c('GO:0000139'),columns = c('SYMBOL'), keytype = "GOALL")
gene_symbols21 <- unique(results$SYMBOL)

a=venn(list(e=genes,b= gene_symbols1	))
b=venn(list(e=genes,b= gene_symbols2	))
c=venn(list(e=genes,b= gene_symbols3	))
d=venn(list(e=genes,b= gene_symbols4	))
e=venn(list(e=genes,b= gene_symbols5	))
f=venn(list(e=genes,b= gene_symbols6	))
g=venn(list(e=genes,b= gene_symbols7	))
h=venn(list(e=genes,b= gene_symbols8	))
i=venn(list(e=genes,b= gene_symbols9	))
j=venn(list(e=genes,b= gene_symbols10	))
k=venn(list(e=genes,b= gene_symbols11	))
l=venn(list(e=genes,b= gene_symbols12	))
m=venn(list(e=genes,b= gene_symbols13	))
n=venn(list(e=genes,b= gene_symbols14	))
o=venn(list(e=genes,b= gene_symbols15	))
p=venn(list(e=genes,b= gene_symbols16	))
q=venn(list(e=genes,b= gene_symbols17	))
r=venn(list(e=genes,b= gene_symbols18	))
s=venn(list(e=genes,b= gene_symbols19	))
t=venn(list(e=genes,b= gene_symbols20	))
u=venn(list(e=genes,b= gene_symbols21	))

print(a)
print(b)
print(c)
print(d)
print(e)
print(f)
print(g)
print(h)
print(i)
print(j)
print(k)
print(l)
print(m)
print(n)
print(o)
print(p)
print(q)
print(r)
print(s)
print(t)
print(u)

gogenes3=read.csv('gocell.csv',header = TRUE)
row.names(gogenes3)=c('MHC class II protein complex',
'MHC protein complex',
'lumenal side of endoplasmic reticulum membrane', 
'integral component of lumenal side of endoplasmic reticulum membrane', 
'ER to Golgi transport vesicle membrane ', 
'coated vesicle membrane', 
'transport vesicle membrane', 
'lytic vacuole membrane', 
'clathrin-coated endocytic vesicle membrane', 
'COPII-coated ER to Golgi transport vesicle',  
'lysosomal membrane',  
'clathrin-coated endocytic vesicle',  
'clathrin-coated vesicle membrane',  
'trans-Golgi network membrane',  
'lysosome',  
'integral component of endoplasmic reticulum membrane',  
'endocytic vesicle membrane',  
'endocytic vesicle', 
'trans-Golgi network', 
'cytoplasmic vesicle membrane', 
'Golgi membrane')

gogenes3[is.na(gogenes3)] <- 0

#heatmap, save image
png("gobio.png", width = 6000, height = 5000,res = 600)
heatmap.2(t(gogenes1), col = blues9,cexCol=0.7,cexRow=0.8,srtCol=7,key=FALSE,trace='none')
dev.off()

png("gomol.png", width = 6500, height = 4000,res = 600)
heatmap.2(t(gogenes2), col = blues9,cexCol=0.6,cexRow=0.8,srtCol=15,key=FALSE,trace='none')
dev.off()

png("gocell.png", width = 6500, height = 4000,res = 600)
heatmap.2(t(gogenes3), col = blues9,cexCol=0.6,cexRow=0.8,srtCol=15,key=FALSE,trace='none')
dev.off()

#Reactome
BiocManager::install("ReactomeContentService4R")

library(ReactomeContentService4R)
R1=event2Ids(event.id = "R-HSA-202430")$geneSymbol
R2=event2Ids(event.id = "R-HSA-202427")$geneSymbol
R3=event2Ids(event.id = "R-HSA-389948")$geneSymbol
R4=event2Ids(event.id = "R-HSA-202433")$geneSymbol
R5=event2Ids(event.id = "R-HSA-388841")$geneSymbol
R6=event2Ids(event.id = "R-HSA-877300")$geneSymbol
R7=event2Ids(event.id = "R-HSA-202424")$geneSymbol
R8=event2Ids(event.id = "R-HSA-2132295")$geneSymbol
R9=event2Ids(event.id = "R-HSA-202403")$geneSymbol
R10=event2Ids(event.id = "R-HSA-913531")$geneSymbol

a=venn(list(e=genes,b= R1	))
b=venn(list(e=genes,b= R2	))
c=venn(list(e=genes,b= R3	))
d=venn(list(e=genes,b= R4	))
e=venn(list(e=genes,b= R5	))
f=venn(list(e=genes,b= R6	))
g=venn(list(e=genes,b= R7	))
h=venn(list(e=genes,b= R8	))
i=venn(list(e=genes,b= R9	))
j=venn(list(e=genes,b= R10	))

print(a)
print(b)
print(c)
print(d)
print(e)
print(f)
print(g)
print(h)
print(i)
print(j)

re=read.csv('re.csv',header = TRUE)
row.names(re)=c('Translocation of ZAP-70 to Immunological synapse',
                'Phosphorylation of CD3 and TCR zeta chains',
                'PD-1 signaling',
                'Generation of second messenger molecules',
                'Costimulation by the CD28 family', 
                'Interferon gamma signaling', 
                'Downstream TCR signaling', 
                'MHC class II antigen presentation',
                'TCR signaling',
                'Interferon Signaling')
                
png("re.png", width = 6500, height = 4000,res = 600)
heatmap.2(t(re), col = blues9,cexCol=0.6,cexRow=0.8,srtCol=15,key=FALSE,trace='none')
dev.off()

library(igraph)
library(RCy3)
cytoscapePing()
renetwork=graph_from_incidence_matrix(re)
createNetworkFromIgraph(renetwork)

go1network=graph_from_incidence_matrix(gogenes1)
createNetworkFromIgraph(go1network)

go2network=graph_from_incidence_matrix(gogenes2)
createNetworkFromIgraph(go2network)

go3network=graph_from_incidence_matrix(gogenes3)
createNetworkFromIgraph(go3network)

#Fisher's Exact Test
BiocManager::install("GeneOverlap")
library(GeneOverlap)

overlapobject=newGeneOverlap(depgenes3[[1]],Agenes3[[1]])
testGeneOverlap(overlapobject)

#Result
#Overlapping p-value=1.5e-05  
