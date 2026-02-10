suppressMessages(library(emmeans))
suppressMessages(library(lme4))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

work.dir <- 'final-data/04.alpha-diveristy/post-hoc-fits'
md <- 
    read.csv(paste(work.dir, 'alpha_trace.tsv', sep='/')
            ,sep='\t'
            ,colClasses='character') %>%
    rename(all_of(c(sample.id = 'X'))) %>%
    mutate(., days_inv = as.numeric(days_inv)
            , observed_features = as.numeric(observed_features)
            , shannon = as.numeric(shannon)
            , simpson = as.numeric(simpson)
            , abx_code2 = as.factor(abx_code2)
            , sex = as.factor(sex)
            , after_travel_na_eu = as.factor(after_travel_na_eu)
            , extraction = as.factor(extraction_no)
            , extraction_no = as.factor(first_extraction)
            , period = as.factor(period)
            , depth_d = as.numeric(depth_d)
            , obs0_lo = as.factor((observed_features_0_cat %in% c(0, '0')) * 1)
            , shannon0_lo = (shannon_0_cat %in% c('0', 0)) * 1
            , simpson0_lo = (simpson_0_cat %in% c('0', 0)) * 1
            ) %>%
    select(-c(contains('0_cat'), observed_features_0, shannon_0, simpson_0))
# Sets up the models
ref.eqs <- list(model0 = 'depth_d'
               ,model1 = 'depth_d + C(extraction_no)'
               ,model2 = 'depth_d + C(extraction_no) + after_travel_na_eu'
               ,model3 = paste('depth_d', 'C(extraction_no)', 
                              'after_travel_na_eu', 'age', 'sex', 
                              'vegetarian', sep=' + ')
               )
metrics <- list(observed_features='observed_features'
               ,shannon='shannon'
               ,simpson='simpson'
               )
random <- list(re = '(1 | host_subject_id)'
              ,rs = '(1 + days_inv | host_subject_id)'
              )
formula <- list()
for(model in names(ref.eqs)) {
    for(met in names(metrics)) {
        for(rs in names(random)) {
            formula[[paste(met, rs, model, sep='.')]] <- 
                paste(metrics[[met]], '~', 'C(abx_code2)*days_inv +',
                      ref.eqs[[model]], '+', random[[rs]], sep=' ')
            }
    }
}

### USEFUL FUNCTIONS ###
# Fits the models using lmer
wrap.lmer <- function(formula){
    fit <- lmer(as.formula(formula), data=md)
    return(fit)
}
# We're going to extract the means from the categorical data and the trends
# for the continous data.
wrap.emmeans.abx <- function(fit.1){
    em <- emmeans(fit.1
                  , pairwise ~  abx_code2 | days_inv
                  , mode = 'satterthwaite'
                  )
    return(em)
}
wrap.etrends.abx <- function(fit.1){
    et <- emtrends(fit.1
                   , pairwise ~ abx_code2
                   , var='days_inv'
                   , mode = 'satterthwaite'
                   )
    return(et)
}
wrap.means.est <- function(em){
    return(as.data.table(em$emmeans))
}
wrap.trends.est <- function(em){
    return(as.data.table(em$emtrends))
}
wrap.means.contrasts <- function(em){
    ctst <- em$contrasts %>%
        as.data.table() %>%
        mutate(., abx0 = str_split_i(contrast, ' - ', 1)
                , abx1 = str_split_i(contrast, ' - ', 2)) %>%
        mutate_at(., .vars=vars(contains('abx'))
                   , .funs=~str_remove(., 'abx_code2')) %>%
    relocate(., abx0, abx1, .after=contrast) %>%
    select(-c(contrast))
    return(ctst)
}

# RUNNING THE DATA
fits <- lapply(formula, wrap.lmer)

# Gets intercept estimate
em.abx <- lapply(fits, wrap.emmeans.abx)
em.abx.beta <-
    lapply(em.abx, wrap.means.est) %>% 
    rbindlist(., id='test') %>%
    mutate(., metric = str_split_i(test, '\\.', 1)
            , rs = str_split_i(test, '\\.', 2)
            , model = str_split_i(test, '\\.', 3)
            , kind = 'b'
            ) %>%
    relocate(metric, model, rs, kind, .before=test) %>%
    rename(all_of(c(beta='emmean'))) %>%
    select(-c(test, days_inv))

em.abx.ctst <- 
    lapply(em.abx, wrap.means.contrasts) %>%
    rbindlist(., id='test') %>%
    mutate(., metric = str_split_i(test, '\\.', 1)
            , rs = str_split_i(test, '\\.', 2)
            , model = str_split_i(test, '\\.', 3)
            , kind = 'b'
            ) %>%
    relocate(metric, model, rs, kind, .before=test) %>%
    select(-c(test, days_inv))


# Gets the slope estimates
et.abx <- lapply(fits, wrap.etrends.abx)
et.abx.beta <- 
    lapply(et.abx, wrap.trends.est) %>% 
    rbindlist(., id='test') %>%
    mutate(., metric = str_split_i(test, '\\.', 1)
            , rs = str_split_i(test, '\\.', 2)
            , model = str_split_i(test, '\\.', 3)
            , kind = 'm'
            ) %>%
    rename(all_of(c(beta='days_inv.trend'))) %>%
    relocate(metric, model, rs, kind, .before=test) %>%
    select(-c(test))
et.abx.ctst <- 
    lapply(et.abx, wrap.means.contrasts) %>%
    rbindlist(., id='test') %>%
    mutate(., metric = str_split_i(test, '\\.', 1)
            , rs = str_split_i(test, '\\.', 2)
            , model = str_split_i(test, '\\.', 3)
            , kind = 'm'
            ) %>%
    relocate(metric, model, rs, kind, .before=test) %>%
    select(-c(test))

# saves the output files
rbind(em.abx.beta, et.abx.beta) %>%
    write.table(., paste(work.dir,  'alpha_expanded_est.tsv', sep='/')
                 , sep='\t'
                 , quote=FALSE)
rbind(em.abx.ctst, et.abx.ctst) %>%
    write.table(., paste(work.dir,  'alpha_expanded_contrasts.tsv', sep='/')
                 , sep='\t'
                 , quote=FALSE)
