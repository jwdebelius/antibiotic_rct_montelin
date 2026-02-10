suppressMessages(library(vegan))
suppressMessages(library(tidyverse))

meta.fp <- 'final-data/08.baseline_compare/meta_baseline.tsv'
dm.fp <- list(
    jaccard = 'final-data/08.baseline_compare/jaccard_baseline.tsv'
   ,braycurtis = './final-data/08.baseline_compare/braycurtis_baseline.tsv'
   ,unweighted_unifrac = './final-data/08.baseline_compare/unweighted_unifrac_baseline.tsv'
   ,weighted_unifrac = './final-data/08.baseline_compare/weighted_unifrac_baseline.tsv'
   )

meta <- read.table(meta.fp
                  ,sep='\t'
                  ,colClasses='character'
                  ,row.name='sample.id', 
                  ,header=TRUE)
meta$depth <- as.numeric(meta$gtdb_depth)
meta$age <- as.numeric(meta$age)
meta$sex <- as.factor(meta$sex)
meta$extraction_no <- as.factor(meta$extraction_no)
meta$antibiotic <- as.factor(meta$antibiotic)
meta$vegetarian <- as.factor(meta$vegetarian)

read.dm <- function(fp){
    dm <- read.table(fp, row.name='X', header=TRUE, sep='\t')
    colnames(dm) <- row.names(dm)
    return(dm)
}

dms <- lapply(dm.fp, read.dm)


fit.models <- function(dm.df){
    md <- meta[row.names(dm.df), ]
    dm <- dist(dm.df)
    crude.abx <- adonis2(dm ~ antibiotic
                         ,data = md
                         ,by = 'terms'
                         ,nperm=999)
    crude.ext <- adonis2(dm ~ extraction_no
                         ,data = md
                         ,by = 'terms'
                         ,nperm=999)
    crude.age <- adonis2(dm ~ age
                         ,data = md
                         ,by = 'terms'
                         ,nperm=999)
    crude.sex <- adonis2(dm ~ sex
                         ,data = md
                         ,by = 'terms'
                         ,nperm=999)
    crude.diet <- adonis2(dm ~ vegetarian
                         ,data = md
                         ,by = 'terms'
                         ,nperm=999)
    full.abx <- adonis2(dm ~ depth + age + sex + extraction_no + vegetarian + antibiotic
                         ,data = md
                         ,by = 'terms'
                         ,nperm=999)
    full.ext <- adonis2(dm ~ depth + age + sex + antibiotic + vegetarian + extraction_no
                         ,data = md
                         ,by = 'terms'
                         ,nperm=999)
    full.age <- adonis2(dm ~ depth + sex + antibiotic + extraction_no + vegetarian + age
                         ,data = md
                         ,by = 'terms'
                         ,nperm=999)
    full.sex <- adonis2(dm ~ depth + antibiotic + extraction_no + vegetarian + age + sex
                         ,data = md
                         ,by = 'terms'
                         ,nperm=999)
    full.diet <- adonis2(dm ~ depth + antibiotic + extraction_no + age + sex + vegetarian
                         ,data = md
                         ,by = 'terms'
                         ,nperm=999)
    p <- c(crude.abx$`Pr(>F)`[1] 
          ,crude.ext$`Pr(>F)`[1]
          ,crude.age$`Pr(>F)`[1]
          ,crude.sex$`Pr(>F)`[1]
          ,crude.diet$`Pr(>F)`[1]
          ,full.abx$`Pr(>F)`[6]
          ,full.ext$`Pr(>F)`[6]
          ,full.age$`Pr(>F)`[6]
          ,full.sex$`Pr(>F)`[6]
          ,full.diet$`Pr(>F)`[6]
          )
    r2 <- c(crude.abx$R2[1] 
           ,crude.ext$R2[1]
          ,crude.age$R2[1]
          ,crude.sex$R2[1]
          ,crude.diet$R2[1]
          ,full.abx$R2[6]
          ,full.ext$R2[6]
          ,full.age$R2[6]
          ,full.sex$R2[6]
          ,full.diet$R2[6]
          )
    models <- c('crude', 'crude', 'crude', 'crude', 'crude',
                'full', 'full', 'full', 'full', 'full')
    key.var <- c('abx', 'extract', 'age', 'sex', 'diet',
                 'abx', 'extract', 'age', 'sex', 'diet')
    summary <- data.frame(model=models
                         ,key_vars = key.var
                         ,R2 = r2
                         ,p_value = p
                         )
    return(summary)
    }

adonis.res <- lapply(dms, fit.models)
i <- 1
for(metric in names(adonis.res)){
    adonis.res[[metric]][, 'metric'] <- metric
    adonis.res[[metric]][, 'nperm'] <- 999
    if(i == 1){
        res.df <- adonis.res[[metric]]
    } else {
        res.df <- rbind(res.df, adonis.res[[metric]])
    }
    i <- i + 1
}

res.df <- res.df[, c('metric', 'model', 'key_vars', 'R2', 'p_value', 'nperm')]

write.table(res.df
           ,'./final-data/06.baseline_compare/adonis_res.tsv'
            ,sep='\t'
            ,quote=FALSE
            ,row.names=FALSE)