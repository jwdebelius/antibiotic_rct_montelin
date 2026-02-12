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
            , obs_feat0_lo = as.factor(observed_features_0_cat)
            , shannon0_lo = as.factor(shannon_0_cat)
            , simpson0_lo = as.factor(simpson_0_cat)
            ) %>%
    select(-c(contains('0_cat'), observed_features_0, shannon_0, simpson_0,
              observed_features, observed_features_0_cat))

# Sets up formula arguments
fixed <-  'depth_d + C(extraction_no) + after_travel_na_eu'
random <- list(re='(1 | host_subject_id)'
               )
metrics <- list(jaccard='jaccard_sc'
               ,unweighted_unifrac='unweighted_unifrac_sc'
               ,braycurtis='braycurtis_sc'
               ,weighted_unifrac='weighted_unifrac_sc'
               )
formulas <- list()
for(re_name in names(random)){
    for(metric in metrics){
        ref_col <- paste('C(', 'shannon', '0_lo)', sep='')
        inter <- paste(ref_col, 'C(abx_code2)', 'C(period)', sep='*')
        rand <- random[[re_name]]
        form.name <- paste(metric, re_name, sep='.')
        formulas[[form.name]] <-
            paste(metric, '~', inter, '+', fixed, '+', rand)
    }
}
wrap.lmer <- 
    function(formula){
        fit <- lmer(as.formula(formula), data=md)
        return(fit)
    }
wrap.emmeans <- 
    function(fit.1){
        em <-  emmeans(fit.1
                       ,pairwise ~ shannon0_lo | (abx_code2 + period)
                       , mode = 'satterthwaite'
                       )
    }
wrap.emmeans.est <- function(em){
        table <- as.data.table(em$emmeans) %>%
        rename_with(., .cols=ends_with('0_lo'), .fn=~'alpha0_lo')
    return(table)
    }
wrap.emmeans.cst <- function(em){
    csts <- as.data.table(em$contrasts)
    return(csts)
    }

beta.ps.lme <- lapply(formulas, wrap.lmer)
beta.ps.emmeans <- lapply(beta.ps.lme, wrap.emmeans)
beta.ps.est <- lapply(beta.ps.emmeans, wrap.emmeans.est) %>%
    rbindlist(., id='metric') %>%
    mutate(., model='model2')
beta.ps.cst <- lapply(beta.ps.emmeans, wrap.emmeans.cst) %>%
    rbindlist(., id='metric') %>%
    mutate(., model='model2')

beta.ps.est %>% 
    write.table(paste(work.dir, 'beta_strata_est.tsv', sep='/')
               ,sep='\t'
               ,quote=FALSE)
beta.ps.cst %>% 
    write.table(paste(work.dir, 'beta_strata_cst.tsv', sep='/')
               ,sep='\t'
               ,quote=FALSE)