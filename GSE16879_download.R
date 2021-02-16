# Load libraries
library(GEOquery)
library(tidyverse)
library(janitor)

#Download GSE dataset 
list <- getGEO("GSE16879")

#Show list details 
show(list)
class(list)
length(list)
names(list)

#Show data elements
data<-list[[1]]

##Create subsetted dataframe with phenotype data only
pd <- pData(data)
names(pd)

##Extract phenotypes of interest
pd %>%
  select(geo_accession, title, `before or after first infliximab treatment:ch1`, `response to infliximab:ch1`, `disease:ch1`, `tissue:ch1`) ->
pd_select

#Get normalized data
ed<-exprs(list[[1]])

#Get feature data
fd<-fData(data)
names(fd)

##Extract features of interest
fd %>%
  select(ID, 'Sequence Source', 'Gene Title', 'Gene Symbol', 'ENTREZ_GENE_ID', 'RefSeq Transcript ID', `Gene Ontology Biological Process`, `Gene Ontology Cellular Component`, `Gene Ontology Molecular Function`) ->
  fd_select

#transpose the expression data to get sample names as the rows
ed %>%
  t() %>%
  as.data.frame() ->
  ed_transposed

#Assign header to first column
ed_transposed %>%
  rownames_to_column(var="geo_accession") ->
  ed_transposed

#Merge Expression and Phenotype daa
left_join(pd_select, ed_transposed, by= "geo_accession") ->
  AntiTNF_expr_pheno

# Write output to file
write.csv(AntiTNF_expr_pheno, file="AntiTNF_expr_pheno.csv")
write.csv(fd_select, file="Transcript Information.csv")

