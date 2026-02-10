import numpy as np
import pandas as pd
import scipy
import statsmodels.api as sms


def extract_recovery_fits(fit, metric, model, slope=False, ref_value='0'):
    """
    Extracts the parameters for a recovery fit
    """
    params = pd.concat(
        axis=1, 
        objs=[fit.params, fit.bse, fit.tvalues, fit.pvalues, fit.conf_int()]
        )
    params.columns = ['params', 'bse', 't_value', 'p_value', 'ci_lo', 'ci_hi']
    # Drops the non-antibiotic rows
    params.drop(index=[i for i in params.index if _check_var(i)],
                inplace=True)
    # Describes the data
    idx = params.index.to_frame()[0]
    params['abx_code'] = \
        params.index.to_frame()[0].apply(find_abx, ref_value=ref_value)
    params['slope'] = \
        params.index.to_frame()[0].apply(lambda x: 'days_inv' in x)
    params['metric'] = metric
    params['model'] = model
    params['rs_slope'] = slope 
    params['ref_value'] = ref_value

    return params


def extract_period_fits(fit, metric, model, slope=False, ref_value='0', 
                        period=None):
    """
    Extracts the parameters for a categorical interaction with period
    """
    # Pulls out hte parameters
    params = pd.concat(
        axis=1, 
        objs=[fit.params, fit.bse, fit.tvalues, fit.pvalues, fit.conf_int()]
        )
    params.columns = ['params', 'bse', 't_value', 'p_value', 'ci_lo', 'ci_hi']
    
    # Drops non-antibiotic data
    params.drop(index=[v for v in params.index if _check_var(v)],
                inplace=True)

    # Describes the data
    params['abx_code'] = \
        params.index.to_frame()[0].apply(find_abx, ref_value=ref_value)
    if period is None:
        params['period'] = \
            params.index.to_frame()[0].apply(lambda x: 'period' in x) + 1
    else:
        params['period'] = period
    params['metric'] = metric
    params['model'] = model
    params['rs_slope'] = slope 
    params['ref_value'] = ref_value

    return params


def extract_abxbl_fit(fit, metric, model, rs_slope, ref_value='0'):
    """
    Extracts the beta estimates for the expaneded fit
    """
    # Pulls out hte parameters
    params = pd.concat(
        axis=1, 
        objs=[fit.params, fit.bse, fit.tvalues, fit.pvalues, fit.conf_int()]
        )
    params.columns = ['params', 'bse', 't_value', 'p_value', 'ci_lo', 'ci_hi']
    
    # Drops non-antibiotic data
    params.drop(index=[v for v in params.index if _check_var(v)],
                inplace=True)

    # Expands the antibiotic and period data
    params['abx_bl'] = [_convert_abx_bl(v) for v in params.index]
    params['abx'] = np.floor(params['abx_bl'] / 10).astype(int)
    params['period'] = params['abx_bl'] - 10 * params['abx']

    # Describes the model
    params['model'] = model
    params['metric'] = metric
    params['ref_value'] = ref_value
    params['rs_slope'] = rs_slope

    return params


def build_tidy_str(x):
    """
    Builds a tidy string
    """
    if str(x['ref_value']) == str(x['abx_code']):
        return '0 (ref)'
    else:
        return '{params:1.2f} ({ci_lo:1.2f}, {ci_hi:1.2f})'.format(**x)

def tidy_pvals(x):
    if x >= 0.01:
        return f'{x:1.2f}'
    elif x >= 0.001:
        return f'{x:1.3f}'
    elif x < 1e-16:
        return '>1e-16'
    else:
        return f'{x:1.0e}'


def find_abx(x, ref_value='0'):
    """
    Looks for the antibiotic code
    """
    if 'abx_code2' in x:
        return x.split(':')[0].split(')[T.')[1].replace(']', '')
    else:
        return str(ref_value)


def _convert_abx_bl(v):
    """
    Extracts the antibiotic period interaction
    """
    if 'Intercept' in v:
        return 1
    else:
        return int(v.replace(']', '').split('[')[-1].replace('T.', ''))


def _check_var(v):
    """
    Checks variables to exclude
    """
    return (('extraction_no' in v) | 
            ('after_travel' in v) | 
            ('depth' in v) | 
            ('age' in v) | 
            ('sex' in v) | 
            ('veg' in v) |
            ('Group Var' in v) | 
            ('Var' in v) | 
            ('Cov' in v)
            )