suppressWarnings(library(phyloseq))
suppressWarnings(library(tidyverse))
suppressWarnings(library(ANCOMBC))
library("frictionless")

# Loads the baseline files
species.fp <- 'final-data/06.baseline_compare/filtered-species-bl.tsv'
meta.fp <- 'final-data/08.baseline_compare/meta_baseline.tsv'

meta.bl <- read.delim(meta.fp, check.names = FALSE, row.names = 1)
meta.bl <- meta.bl[row.names(meta.bl) != "#q2:types", ]
meta.bl$age <- as.numeric(meta.bl$age)
meta.bl$sex <- as.factor(meta.bl$sex)
meta.bl$abx_code2 <- case_when(meta.bl$abx_code2 %in% c('0', 0) ~ 'A'
                              ,meta.bl$abx_code2 %in% c('1', 1) ~ 'B'
                              ,meta.bl$abx_code2 %in% c('2', 2) ~ 'C'
                              ,meta.bl$abx_code2 %in% c('3', 3) ~ 'D'
                              ,meta.bl$abx_code2 %in% c('4', 4) ~ 'E'
                              ,meta.bl$abx_code2 %in% c('5', 5) ~ 'F'
                              )
meta.bl$abx_code2 <- as.factor(meta.bl$abx_code2)
meta.bl$extraction_no <- as.factor(meta.bl$extraction_no)
meta.bl$antibiotic <- as.factor(meta.bl$antibiotic)
meta.bl$vegetarian <- as.factor(meta.bl$vegetarian)

md.bl <- sample_data(meta.bl)
row.names(md.bl) <- rownames(meta.bl)

otus <- read.delim(species.fp, check.names = FALSE, row.names = 1)
taxa <-  otus[,c('labels')]
otus <- otus[, rownames(meta.bl)]
otu.bl <- otu_table(otus, taxa_are_rows = TRUE)
otu.data.bl <- phyloseq(otu.bl, md.bl)

fit.otus.bl <- ancombc2(data=otu.data.bl,
                        fix_formula = 'age + sex + vegetarian + abx_code2 + extraction_no', 
                        alpha = 0.05)
res.otu.bl <- fit.otus.bl$res %>%
    tibble() %>%
    select(c(taxon, contains('abx_code')))

res.otu.bl %>%
    write.table('final-data/06.baseline_compare/antibiotic-taxa-ancombc-res.tsv'
                ,sep='\t'
                ,quote=FALSE)