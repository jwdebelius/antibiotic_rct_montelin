suppressMessages(library(emmeans))
suppressMessages(library(lme4))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

work.dir <- 'final-data/06.beta-diveristy/post-hoc-fits'
md <- read.csv(paste(work.dir, 'beta_traces.tsv', sep='/')
              ,sep='\t'
              ,colClasses='character') %>%
    data.table() %>%
    mutate_at(., .vars=vars(contains('sc'))
               , .funs=~as.numeric(.) / 10
               ) %>%
    mutate(., days_inv = as.numeric(days_inv)
            , observed_features = as.numeric(observed_features)
            , shannon = as.numeric(shannon)
            , simpson = as.numeric(simpson)
            , age = as.numeric(age)
            , abx_code2 = as.factor(abx_code2)
            , sex = as.factor(sex)
            , after_travel_na_eu = as.factor(after_travel_na_eu)
            , extraction = as.factor(extraction_no)
            , extraction_no = as.factor(first_extraction)
            , period = as.factor(period)
            , depth_d = as.numeric(depth_d)
            )

# Sets up the models
ref.eqs <- list(model0 = 'depth_d'
               ,model1 = 'depth_d + C(extraction_no)'
               ,model2 = 'depth_d + C(extraction_no) + after_travel_na_eu'
               ,model3 = paste('depth_d', 'C(extraction_no)', 
                              'after_travel_na_eu', 'age', 'sex', 
                              'vegetarian', sep=' + ')
                )
metrics <- list(jaccard='jaccard_sc'
               ,unweighted_unifrac='unweighted_unifrac_sc'
               # ,rpca='rpca_sc'
               # ,phylo_rpca='phylo_rpca_sc'
               ,braycurtis='braycurtis_sc'
               ,weighted_unifrac='weighted_unifrac_sc'
               )
random <- list(re = '(1 | host_subject_id)'
              ,rs = '(1 + period | host_subject_id)'
              )
formula <- list()
for(model in names(ref.eqs)) {
    for(met in names(metrics)) {
        for(rs in names(random)) {
            formula[[paste(met, rs, model, sep='.')]] <- 
                paste(metrics[[met]], '~', 'C(abx_code2)*C(period) +',
                      ref.eqs[[model]], '+', random[[rs]], sep=' ')
            }
    }
}

### USEFUL FUNCTIONS ###
wrap.lmer <- 
    function(formula){
        fit <- lmer(as.formula(formula), data=md)
        return(fit)
    }
wrap.emmeans <- 
    function(fit.1){
        em <-  emmeans(fit.1
                       , ~ abx_code2 * period
                       , mode = 'satterthwaite'
                       )
    }
wrap.contrasts <- 
    function(em2){
        csts <- contrast(em2, 'pairwise', simple='each', combine=TRUE)
    }

# Fit data
fits <- lapply(formula, wrap.lmer)
est.m <- lapply(fits, wrap.emmeans)
cnstrst <- lapply(est.m, wrap.contrasts)

# Tidies the estimates4
est.m.2 <- 
    lapply(est.m, as.data.frame) %>%
    rbindlist(., id='test') %>%
    mutate(., metric = str_split_i(test, '\\.', 1)
            , rs = str_split_i(test, '\\.', 2)
            , model = str_split_i(test, '\\.', 3)
            ) %>%
    relocate(metric, model, rs, .before=test) %>%
    rename(all_of(c(beta='emmean'))) %>%
    select(-c(test))

# Tidies contrasts
cst.combo <- 
    lapply(cnstrst, as.data.frame) %>%
    rbindlist(., id='test') %>%
    mutate(., c1 = str_split_i(contrast, " - ", 1)
            , c2 = str_split_i(contrast, ' - ', 2)
            , metric = str_split_i(test, '\\.', 1)
            , rs = str_split_i(test, '\\.', 2)
            , model = str_split_i(test, '\\.', 3)
            ) %>%
    mutate_at(., .vars=vars(c(c1, c2))
               , .funs=~str_remove(., 'abx_code2')
              ) %>%
    mutate_at(., .vars=vars(c(c1, c2))
               , .funs=~str_remove(., 'period')
              ) %>%
    relocate(c1, c2, .after=contrast) %>%
    relocate(metric, model, rs, .before=test) %>%
    select(-c(test, contrast))
    
# Saves output files
write.table(est.m.2
           ,paste(work.dir, 'beta_expanded_est.tsv', sep='/')
           ,sep='\t'
           ,quote=FALSE
           )
write.table(cst.combo
           ,paste(work.dir, 'beta_expanded_contrasts.tsv', sep='/')
           ,sep='\t'
           ,quote=FALSE
           )