suppressMessages(library(emmeans))
suppressMessages(library(lme4))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

# Sets up the model
work.dir <- 'final-data/04.alpha-diveristy/post-hoc-fits'
md <- 
    read.csv(paste(work.dir, 'alpha_trace.tsv', sep='/')
            ,sep='\t'
            ,colClasses='character') %>%
    rename(all_of(c(sample.id = 'X'))) %>%
    mutate(., days_inv = as.numeric(days_inv)
            , obs_feat = as.numeric(observed_features)
            , shannon = as.numeric(shannon)
            , simpson = as.numeric(simpson)
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
               ,rs='(1 + days_inv | host_subject_id)'
               )

metrics <- c('obs_feat'
            ,'shannon'
            ,'simpson')

formulas <- list()
for(re_name in names(random)){
    for(metric in metrics){
        ref_col <- paste('C(', 'shannon', '0_lo)', sep='')
        inter <- paste(ref_col, 'C(abx_code2)', 'days_inv', sep='*')
        rand <- random[[re_name]]
        form.name <- paste(metric, re_name, sep='.')
        formulas[[form.name]] <-
            paste(metric, '~', inter, '+', fixed, '+', rand)
    }
}

fit.model <- function(formula){
    return(lmer(as.formula(formula), data=md))
}
# Fits the model
fits.lme <- lapply(formulas, fit.model)

fit.emmeans.strat <- function(fit){
    em.res <- emmeans(fit
                     ,pairwise ~ shannon0_lo | (abx_code2 + days_inv)
                     ,mode = 'satterthwaite')
    return(em.res)
}
fit.etrends.start <- function(fit){
    em.res <- emtrends(fit
                     ,pairwise ~ shannon0_lo | abx_code2
                     ,var='days_inv'
                     ,mode = 'satterthwaite')
    return(em.res)
}

# Gets the trends for the data
emmeans.strat <- lapply(fits.lme, fit.emmeans.strat)
emmeans.trends <- lapply(fits.lme, fit.etrends.start)
wrap.means.est <- function(em){
    table <- as.data.table(em$emmeans) %>%
        rename_with(., .cols=ends_with('0_lo'), .fn=~'alpha0_lo')
    return(table)
    }
wrap.trends.est <- function(em){
    table <- as.data.table(em$emtrends) %>%
        rename_with(., .cols=ends_with('0_lo'), .fn=~'alpha0_lo')
    return(table)
    }


wrap.means.cst <- function(em){
    csts <- as.data.table(em$contrasts)
    return(csts)
}
emmeans.est <- 
    lapply(emmeans.strat, wrap.means.est) %>%
    rbindlist(., id='metric') %>%
    select(-c(days_inv)) %>%
    rename(all_of(c(`estimate` = 'emmean'))) %>%
    mutate(., model = 'model2'
            #, rs = str_split_i(metric, '.', 2)
            #, metric = str_split_i(metric, '.', 1)
            , kind = 'b'
            )

emmeans.cst <- 
    lapply(emmeans.strat, wrap.means.cst) %>%
    rbindlist(., id='metric') %>%
    select(-c(days_inv, contrast)) %>%
    mutate(., model = 'model2'
            #, rs = str_split_i(metric, '.', 2)
            #, metric = str_split_i(metric, '.', 1)
            , kind = 'b')

etrends.est <- 
    lapply(emmeans.trends, wrap.trends.est) %>%
    rbindlist(., id='metric') %>%
    rename(all_of(c(`estimate` = 'days_inv.trend'))) %>%
    mutate(., model = 'model2'
            #, rs = str_split_i(metric, '.', 2)
            #, metric = str_split_i(metric, '.', 1)
            , kind = 'm')

etrends.csts <- 
    lapply(emmeans.trends, wrap.means.cst) %>%
    rbindlist(., id='metric') %>%
    select(-c(contrast)) %>%
    mutate(., model = 'model2'
            #, rs = str_split_i(metric, '.', 2)
            #, metric = str_split_i(metric, '.', 1)
            , kind = 'm')


rbind(emmeans.est, etrends.est) %>%
    write.table(., paste(work.dir,  'alpha_strat_est.tsv', sep='/')
                 , sep='\t'
                 , quote=FALSE)

rbind(emmeans.cst, etrends.csts)%>%
    write.table(., paste(work.dir,  'alpha_strat_contrasts.tsv', sep='/')
                 , sep='\t'
                 , quote=FALSE)
