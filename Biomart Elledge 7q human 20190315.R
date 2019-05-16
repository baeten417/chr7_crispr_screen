setwd("~/data/Protocols/Baeten_data/data mining")
a <- read.delim("~/data/Protocols/Baeten_data/data mining/Elledge 2018.txt")
library("biomaRt")
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters=listFilters(ensembl)
head(filters)
attributes=listAttributes(ensembl)
head(attributes)

b <- getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position"), mart=ensembl)
head(b)
write.csv(b, file="gene_chr_start.out")

c <- read.delim("~/data/analysis/reference/hg19/gene_chr_start_hg38.txt")
a <- read.delim("~/data/Protocols/Baeten_data/data mining/Elledge 2018.txt")
m <- merge(a, c, by="gene")
head(m)
write.csv(m, file="Elledge_chr.out")

###Make AML and HSC reference gene sets
a <- read.delim("AMLexpressedGenes.txt")
hist(a$del.7q..FPKM, breaks=100)
hist(a$otherAML..FPKM, breaks=100)
length(which(a$del.7q..FPKM > 10))
length(which(a$otherAML..FPKM > 10))

a <- read.delim("Elledge_chr7_forR.txt")
aml <- read.delim("AMLexpressedGenes.txt")
hsc <- read.delim("HSC_RNA-seq.txt")
m <- merge(a, aml, by="gene", all.x=T)
m <- merge(m, hsc, by="gene", all.x=T)
write.csv(m, file="Elledge_chr7_HSC_AML.out")

###Get gene descriptions
searchAttributes(mart=ensembl, pattern="description")
searchFilters(mart=ensembl, pattern="description")

out <- getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "description"), mart=ensembl)

head(out)
write.csv(out, file="gene_coordinate_description.out")

IDs <- c("CUX1", "NRF1")
genedesc <- getBM(attributes=c('external_gene_name','description'), filters = 'external_gene_name', values = IDs, mart =ensembl)
genedesc
test <- getBM(attributes='wikigene_description', filters="external_gene_name", values=IDs, mart=ensembl)
test
