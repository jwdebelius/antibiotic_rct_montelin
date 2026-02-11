suppressWarnings(library(phyloseq))
suppressWarnings(library(tidyverse))
suppressWarnings(library(ANCOMBC))
# library("optparse")
# library("frictionless")
# library("jsonlite")

# Sets up the filepaths necessary
table.in.fp <- 'final-data/07.ancom-bc/table_for_ancom.tsv'
meta.in.fp <- 'final-data/01.metadata/analysis_sample_meta.tsv'

# ### md ###
# Sets up the md file and codes the key data as factors
md <- read.delim(meta.in.fp, check.names = FALSE, row.names = 1)
md <- md[row.names(md) != "#q2:types", ]
# We code the numeric data as numbers, and then the two covariates I care
# about as factors
md$days_since_intervention <- as.numeric(md$days_since_intervention)
md$days_inv <- 30 / md$days_since_intervention
md$extraction_no <- as.factor(md$extraction_no)
md$abx_code2 <- as.factor(md$abx_code2)
#md$host_subject_id <- as.factor(md$host_subject_id)
md$after_travel_na_eu <- as.factor(md$after_travel_na_eu)

# Creates dummy interaction because ANCOM-BC keeps breaking in new and 
# intersting values
md$days_inv_1 = (md$days_inv)*(md$abx_code2 == '1')
md$days_inv_2 = (md$days_inv)*(md$abx_code2 == '2')
md$days_inv_3 = (md$days_inv)*(md$abx_code2 == '3')
md$days_inv_4 = (md$days_inv)*(md$abx_code2 == '4')
md$days_inv_5 = (md$days_inv)*(md$abx_code2 == '5')

meta.fu <- md[md$days_since_intervention > 0, ]
# Table
table <- read.delim(table.in.fp
                    , check.names = FALSE
                    , row.names = 1)
table.fu <- table[,row.names(meta.fu)]

# Follow up curves
md.fu <- sample_data(meta.fu)
row.names(md.fu) <- rownames(meta.fu)
otu.fu <- otu_table(table.fu, taxa_are_rows=TRUE)
data.fu <- phyloseq(otu.fu, md.fu)

fit.tax.fu <- ancombc2(
    data = data.fu, 
    fix_formula = 'abx_code2 + extraction_no + after_travel_na_eu + days_inv + days_inv_1 + days_inv_2 + days_inv_3 + days_inv_4 + days_inv_5', 
    rand_formula = '(1 | host_subject_id)',
    alpha = 0.05,
    lme_control = lme4::lmerControl(),
    )
print('res done')
write.table(fit.tax.fu$res, 'final-data/07.ancom-bc/follow_up_res.tsv',
            sep='\t', quote=FALSE)
