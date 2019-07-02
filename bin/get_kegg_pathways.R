#! /usr/bin/env Rscript

# Get KEGG pathways from hsa identifiers.

###################
# LIBRARIES
###################

library(KEGGREST)
res <- keggLink("pathway", "hsa")
write.table(cbind(names(res), res), col.names=FALSE, quote=FALSE, sep="\t", row.names=FALSE, file = "kegg_pathways.txt")

