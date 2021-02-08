#! /usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

file <- args[1]
column <- as.integer(args[2])
output <- args[3]

data <- read.table(file, header = FALSE , sep="\t")

pdf(output)

obj <- ggplot(data, aes(x=data[[column]]))
obj <- obj + geom_density()
obj

dev.off()
