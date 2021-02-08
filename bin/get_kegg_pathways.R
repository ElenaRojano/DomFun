#! /usr/bin/env Rscript

# Get KEGG pathways from hsa identifiers.

args <- commandArgs(trailingOnly = TRUE)
kegg_organism <- args[1]
output_file <- args[2]

###################
# LIBRARIES
###################

library(KEGGREST)
res <- keggLink("pathway", kegg_organism)
write.table(cbind(names(res), res), col.names=FALSE, quote=FALSE, sep="\t", row.names=FALSE, file = output_file)