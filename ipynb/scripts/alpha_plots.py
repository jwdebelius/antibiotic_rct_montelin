import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import ticker
import seaborn as sn

abx_order = ['untreated', 'nitrofurantoin', 'pivmecillinam',
             'co-trimoxazole', 'ceftibuten', 'ciprofloxacin']
colors = ['#969696', '#009E73','#9D654C', 
          '#4E7705', '#098BD9',  '#CC79A7']
spectrum_order = ['untreated', 'narrow', 'broad']
alpha_metrics = ['observed_features', 'shannon', 'simpson']
abx_names = {str(i): abx for i, abx in enumerate(abx_order)}
abx_colors = {str(i): color for i, color in enumerate(colors)}

clean_metrics = {'observed_features': 'Observed Species',
                 'unweighted_unifrac': 'unweighted\nUniFrac',
                 'weighted_unifrac': 'weighted\nUniFrac',
                 'braycurtis': 'Bray Curtis',
                 'rpca': 'rPCA\nAitchison',
                 'phylo_rpca': 'phylo rPCA\nAitchison'
                 }
abx_colors2 = {i: color for i, color in enumerate(colors)}

hide_axes = dict(left=False, labelleft=False,
                 right=False, labelright=False,
                 length=0,
                 labelsize=0)


def _add_labelfig(figya, label, fontsize=10, rotation=0, 
                  position=(0.05, 0.1, 0.9, 0.85)):
    """
    Creates a figure panel wtih the xlabel
    """
    axy = figya.add_subplot(1,1,1, 
                            position=position, 
                            facecolor='None')
    axy.xaxis.set_tick_params(**hide_axes)
    axy.yaxis.set_tick_params(**hide_axes)
    axy.text(0.5, 0.5, 
             label, 
             size=fontsize, 
             rotation=rotation, 
             ha='center', 
             va='center')
    sn.despine(figya, left=True, right=True, top=True, bottom=True)


def _plot_forest_square(ax: plt.Axes, 
                        abx: str, 
                        sum_df: pd.DataFrame, 
                        xx : np.array=np.arange(0, 7),
                        error_kws: dict=None, 
                        shade_kws: dict=None):
    """
    Makes a forest plot square
    """
    error_kwargs = dict(capsize=3,
                        marker='.',
                        linestyle='',
                        zorder=3,
                        )
    if isinstance(error_kws, dict):
        error_kwargs.update(error_kws)

    shade_kwargs = dict(alpha=0.25,
                        linewidth=0,
                        zorder=0,
                        )
    if isinstance(shade_kws, dict):
        shade_kwargs.update(shade_kws)
    
    if abx > 0:
        ax.fill_between(xx,
                        sum_df.loc[0, 'mean'] - sum_df.loc[0, 'ci'],
                        sum_df.loc[0, 'mean'] + sum_df.loc[0, 'ci'],
                        color=colors[0],
                        **shade_kwargs,
                        )
    ax.errorbar(x=xx,
                y=sum_df.loc[abx, 'mean'],  
                yerr=sum_df.loc[abx, 'ci'],
                color=colors[int(abx)],
                **error_kwargs,
                )


def _plot_trace_square(ax: plt.Axes, 
                       abx: str, 
                       sum_df: pd.DataFrame, 
                       trace_df: pd.DataFrame, 
                       xx : np.array=np.arange(0, 7),
                       shade_kws=None, 
                       mean_kws=None, 
                       trace_kws=None, 
                       scale=0.05,
                       ):
    """
    Plots a per individual trace
    """
    # Shade kws
    shade_kwargs = dict(color=colors[abx],
                        alpha=0.25,
                        linewidth=0,
                        zorder=0,)
    if isinstance(shade_kws, dict):
        shade_kwargs.update(shade_kws)
    # Mean trace kwargs
    mean_kwargs = dict(linestyle='-',
                       marker='.',
                       color=colors[abx],
                       zorder=2,)
    if isinstance(mean_kws, dict):
        mean_kwargs.update(mean_kws)
    # Trace_kwargs 
    trace_kwargs = dict(linestyle='-',
                        marker='.',
                        linewidth=1/3,
                        markersize=1,
                        color='#525252')
    if isinstance(trace_kws, dict):
         trace_kwargs.update(trace_kws)

    # Gloabl alpha diveristy
    ymean = sum_df.loc[abx, 'mean']
    yerr = sum_df.loc[abx, 'ci']
    ax.fill_between(x=xx, y1=ymean - yerr, y2=ymean + yerr,
                    **shade_kwargs) 
    ax.plot(xx, ymean, **mean_kwargs)
    

    trace_x = np.vstack([trace_df.columns.astype(float)] * len(trace_df)) + \
        np.vstack([np.random.normal(size=len(trace_df))] * len(xx)).T * scale
    
    ax.plot(trace_x.T, trace_df.T, **trace_kwargs)


def _format_alpha_axes(ax, abx, bottom):
    """
    Lazy function to clean up axes
    """
    ax.set_xticks(np.arange(7))
    ax.set_xlim(-0.5, 6.5)
    ax.set_xticklabels([-1, 7, 14, 30, 90, 180, 365], size=9, 
                       rotation=90, ha='center', va='top')
    ax.xaxis.set_tick_params(bottom=False, right=False, length=0,
                             labelbottom=bottom)
    ax.yaxis.set_tick_params(left=True, labelleft=(abx==0),
                             right=False, labelsize=9)
    ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:1.0f}%'))


def _format_beta_axes(ax, abx, bottom):
    """
    Lazy function to clean up axes
    """
    ax.set_xticks(np.arange(6) + 1)
    ax.set_xlim(0.5, 6.5)
    ax.set_xticklabels([7, 14, 30, 90, 180, 365], size=9, 
                       rotation=90, ha='center', va='top')
    ax.xaxis.set_tick_params(bottom=False, right=False, length=0,
                             labelbottom=bottom)
    ax.yaxis.set_tick_params(left=True, labelleft=(abx==0),
                             right=False, labelsize=9)


def _plot_betas(axf, plot_df, i, j, error_kws=None):
    """
    Builds a forest plot with an inset
    """
    error_kwargs = dict(linestyle='',
                        capsize=3,
                        markersize=5,
                        )
    if isinstance(error_kws, dict):
        error_kwargs.update(error_kws)

    # adds the axis labels as appropriate
    axf.plot([0, 0], [-1, 6], 'k-', linewidth=1, zorder=0)
    axf.set_ylim(-0.5, 5.5)
    axf.set_yticks(plot_df['y'])
    axf.set_yticklabels(plot_df['tidy'], ha='center')
    axf.set_yticks(5 - np.arange(6), minor=True)
    axf.set_yticklabels(abx_order, minor=True)

    # Forest plot
    for (mark, c, mfc), df_ in plot_df.groupby(['marker', 'color', 'facecolor']):
        axf.errorbar(x=df_['params'],
                     xerr=df_['width'],
                     y=df_['y'],
                     color=c,
                     marker=mark,
                     markerfacecolor=mfc,
                     **error_kwargs
                     )
    sn.despine(ax=axf, left=True, right=True, top=True, bottom=False)
    axf.text(axf.get_xlim()[0], axf.get_ylim()[1], 'ABCDEFG'[int(3 * j + i)],
             ha='left', va='top', size=12)

    axf.xaxis.set_tick_params(bottom=True, labelbottom=True, labelsize=8,
                              top=False, labeltop=False)
    axf.yaxis.set_tick_params(which='major', left=False, labelleft=False, 
                              right=False, labelright=True, labelsize=8, 
                              length=30)
    axf.yaxis.set_tick_params(which='minor', left=False, labelleft=True,
                              right=False, labelright=False, labelsize=8,
                              length=0)
    sn.despine(ax=axf, left=True, right=True, top=True, bottom=False)


def plot_caption_with_model(fig, color='#252525'):
    """
    Helper function to make a legend for model 2 and model 3
    """
    # Marker kwargs
    marker_kws = dict(xerr = 0.025, capsize=3, color=color, markersize=5)
    
    ax = fig.add_subplot(1,1,1, facecolor='None')
    ax.xaxis.set_tick_params(**hide_axes)
    ax.yaxis.set_tick_params(**hide_axes)
    ax.set_ylim(-0.5, 2.5)
    ax.set_xlim(0, 1)
    
    # Error bar
    ax.errorbar(x = [0.05], y=2, marker='s', **marker_kws)
    ax.errorbar(x = [0.05], y=1, marker='o', **marker_kws)
    ax.errorbar(x = [0.05], y=0, marker='v', markerfacecolor='w', **marker_kws)
    ax.errorbar(x = [0.30], y=0, marker='v', **marker_kws)

    # Labes
    ax.text(x=0.1, y=2, ha='left', va='center', size=8,
            s=('Model 1: Adjusted for difference in depth, extraction plate, '
               'and after travel outside North America or the EU'))
    ax.text(x=0.1, y=1, ha='left', va='center', size=8,
            s=('Model 2: Model 1 plus adjustment for age, sex, and '
               'vegetarian diet'))
    ax.text(x=0.1, y=0, ha='left', va='center', size=8, s=('p ≥ 0.05'))
    ax.text(x=0.35, y=0, ha='left', va='center', size=8, s=('p < 0.05'))



def plot_caption_with_slope(fig, color='#252525'):
    """
    Helper function to make a legend for model 2, random intercept/slope
    """
    # Marker kwargs
    marker_kws = dict(xerr = 0.025, capsize=3, color=color, markersize=5)
    
    ax = fig.add_subplot(1,1,1, facecolor='None')
    ax.xaxis.set_tick_params(**hide_axes)
    ax.yaxis.set_tick_params(**hide_axes)
    ax.set_ylim(-0.5, 2.5)
    ax.set_xlim(0, 1)

    # Markers
    ax.errorbar(x = [0.05], y=1, marker='s', **marker_kws)
    ax.errorbar(x = [0.50], y=1, marker='D', **marker_kws)
    ax.errorbar(x = [0.05], y=0, marker='v', markerfacecolor='w', **marker_kws)
    ax.errorbar(x = [0.50], y=0, marker='v', **marker_kws)
    
    # Text
    model_ = ('Model adjusted for difference in depth from baseline, '
              'extraction plate, and after travel outside North America or '
              'the EU')
    ax.text(x=0, y=2, s=model_, ha='left', va='center', size=8)
    ax.text(x=0.10, y=1, s='Random intercept', ha='left', va='center', size=8)
    ax.text(x=0.55, y=1, s='Random intercept and slope', 
            ha='left', va='center', size=8)
    ax.text(x=0.1, y=0, ha='left', va='center', size=8, s=('p ≥ 0.05'))
    ax.text(x=0.55, y=0, ha='left', va='center', size=8, s=('p < 0.05'))


def _expand_cols(params, timing, col):
    """Converts the expanded data to a timing curve"""
    exp_df = params[[col] * len(timing.columns)].values * timing
    exp_df.reset_index(inplace=True)
    exp_df = exp_df.melt(id_vars=['metric', 'model', 'rs', 'slope', 'abx_code2'],
                         var_name='days',
                         value_name=col
                         )         
    exp_df = \
        exp_df.groupby(['metric', 'model', 'rs', 'abx_code2', 'days'])[col].sum()
    return exp_df


def build_timing_grid(days, params):
    """
    Calculates the timing above and below for the data set

    This would probably be better calculated with the var-covar matrix
    but this is probably a good enough approximation
    """

    timing = pd.DataFrame(data=np.vstack([30 / days] * len(params)),
                          columns=days,
                          index=params.index)
    # Covers the intercept
    timing.mask(timing.index.to_frame()['slope'].isin([0, '0']), 1, 
                inplace=True)

    exp_ = pd.DataFrame({col_: _expand_cols(params, timing, col_)
                         for col_ in ['beta', 'width']})
    exp_['timing'] = exp_.index.to_frame()['days']
    exp_['ci_lo'] = exp_['beta'] - exp_['width']
    exp_['ci_hi'] = exp_['beta'] + exp_['width']

    return exp_


def _plot_trace_ax(ax, abx, data, dash_thresh=7, line_kws=None):
    """
    Plots the antibiotic recovery curve
    """
    line_kwargs = dict(linewidth=1.5)
    if isinstance(line_kws, dict):
        line_kwargs.update(line_kws)
    line_kwargs = {k:v for k, v in line_kwargs.items()
                   if k not in {'color', 'linestyle'}}
    
    df = data.xs(abx, level='abx_code2')

    dashed = df.loc[(df['timing'] <= dash_thresh)].index
    solid = df.loc[(df['timing'] >= dash_thresh)].index

    ax.fill_between(x='timing', y1='ci_lo', y2='ci_hi', data=df,
                    alpha=0.25,
                    color=abx_colors[abx],
                    edgecolor='None')
    ax.plot(df.loc[dashed, 'timing'], 
            df.loc[dashed, 'beta'],
            color=abx_colors[abx], 
            linestyle=':',
            **line_kwargs
            )
    ax.plot(df.loc[solid, 'timing'], 
            df.loc[solid, 'beta'],
            color=abx_colors[abx], 
            linestyle='-',
            **line_kwargs
            )


def _build_beta_boxplot(ax, metric, abx, meta_r):
    """
    Builds the paired boxplot and beta fits
    """
    box_kwargs = dict(
        data=meta_r.loc[meta_r['abx_code2'].isin([abx, 0])],
        x='period',
        order=[1, 2],
        ax=ax,
        linewidth=1,
        flierprops=dict(marker='o', markerfacecolor='None', markersize=3, 
                        markeredgewidth=0.5),
        )
    if abx in {'0', 0}:
        sn.boxplot(y=metric,                    
                   color=colors[0],
                   **box_kwargs,
                   )
    else:
        sn.boxplot(y=metric, 
                   hue='abx_code2',
                   hue_order=[0, abx],
                   palette=abx_colors2,
                   **box_kwargs,
                   width={0: 0.6, '0': 0.6}.get(abx, 0.8)
                   )
        ax.get_legend().set_visible(False)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['7-30', '≥90'])
    ax.set_ylabel(clean_metrics.get(metric, metric.title()), size=9)
    ax.yaxis.set_tick_params(left=True, labelleft=(abx == 0),
                             right=False, labelright=False,
                             labelsize=8)
    
    ax.xaxis.set_tick_params(bottom=True, labelbottom=True, 
                             top=False, labeltop=False,
                             labelsize=8)
    ax.xaxis.get_label().set_visible(False)
    ax.yaxis.get_label().set_visible(abx == 0)


def _sub_upper(y1, y2, ydiff, ypad=1, top_pad=1):
    """
    Puts together the line spacing
    """
    upper = np.max([y1, y2]) + top_pad * ydiff
    tops = np.array([y1, upper, upper, y2]) + ypad * ydiff

    return tops


def _build_lines_abx(ylim, abx, metric, tops, period_pvals, inter_pvals,
                    xscale=0.2, yscale=0.05, ypad=1, 
                     top_pad=1.5):
    """
    Generates the plotting lines for marking significance
    """
    ydiff = np.diff(ylim)[0] * yscale
    if abx in {'0', 0}:
        y_vs = _sub_upper(*tops.loc[[(int(abx), 1), (int(abx), 2)], metric], 
                          ydiff, ypad, top_pad)
        x_vs = np.array([0, 0, 1, 1])
        p_str = inter_pvals.loc[(metric, '0')]
        return ([[x_vs, y_vs]], [p_str])
    else:
        y_p1 = _sub_upper(*tops.loc[[(0, 1), (abx, 1)], metric],
                          ydiff, ypad, top_pad)
        y_p2 = _sub_upper(*tops.loc[[(0, 2), (abx, 2)], metric],
                          ydiff, ypad, top_pad)
        y_vs = _sub_upper(y_p1.max(), y_p2.max(), ydiff, ypad * 2, top_pad)

        x_p1 = np.array([-1, -1, 1, 1]) * xscale + 0
        x_p2 = np.array([-1, -1, 1, 1]) * xscale + 1
        x_vs = np.array([0, 0, 1, 1])

        p_str = np.hstack([
            period_pvals.loc[[(str(abx), 1), (str(abx), 2)], metric], 
            inter_pvals.loc[(metric, str(abx))]])

        return [[x_p1, y_p1], [x_p2, y_p2], [x_vs, y_vs]] , p_str


def _wrap_striped_forest(ax, metric, metric_df, labelleft=False):
    """
    Wraps the forest plot clean up
    """
    # Builds the forest plot
    for (c, m, fc), df_ in metric_df.groupby(['color', 'marker', 'fc']):
        ax.errorbar(y=df_['y2'],
                    x=df_['params'],
                    xerr=df_['width'],
                    linestyle='',
                    marker=m,
                    color=c,
                    markerfacecolor=fc,
                    capsize=3,
                    markersize=4,
                    markeredgewidth=1,
                   )
    # Plots the reference line and timing reference
    ax.plot([0, 0], [-2, 15], 'k-', linewidth=1, zorder=1)
    ax.set_ylim(-1, 12.5)
    xlim = ax.get_xlim()
    for y_ in metric_df.groupby('abx_code')['y2'].mean():
        ax.fill_between(x=np.array(xlim) * 1.5, 
                        y1=y_, y2=y_-1, 
                        color='#f1f1f1', 
                        zorder=0)
    ax.set_xlim(xlim)
    
    # Sets up the antibiotic label on the yaxis
    ax.set_yticks(metric_df.groupby('abx_code')['y2'].mean(), minor=False)
    ax.set_yticklabels(abx_order, minor=False)
    ax.yaxis.set_tick_params(which='major',
                             left=False, labelleft=labelleft,
                             right=False, labelright=False,
                             labelsize=8, length=0)
    
    # adds minor label as tidy
    ax.set_yticks(metric_df['y2'], minor=True)
    ax.set_yticklabels(metric_df['tidy'], 
                       minor=True,
                       ha='center', 
                       weight='light')
    ax.yaxis.set_tick_params(which='minor',
                             left=False, labelleft=False,
                             right=False, labelright=True,
                             length=32, labelsize=8,
                             )

    # Cleans up the x-axis
    ax.xaxis.set_tick_params(bottom=True, labelbottom=True,
                         top=False, labeltop=False,
                         labelsize=8)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(4))
    ax.set_xlabel('$\\Delta_{baseline}$', size=10)

    sn.despine(ax=ax, left=True, right=True, top=True, bottom=False)


def _wrap_striped_legend_model(figl):
    """
    Wraps the forest legend
    """
    gsl = figl.add_gridspec(1, 2, width_ratios=[8, 1])
    axl1 = figl.add_subplot(gsl[0], facecolor='None')
    axl1.errorbar(x = [0.05], y=2, xerr=0.025, marker='s', capsize=3, color='#252525', markersize=5)
    axl1.errorbar(x = [0.05], y=0.8, xerr=0.025, marker='o', capsize=3, color='#252525', markersize=5)
    axl1.errorbar(x = [0.05], y=0, xerr=0.025, marker='v', capsize=3, color='#252525', markerfacecolor='w', markersize=5)
    axl1.errorbar(x = [0.30], y=0, xerr=0.025, marker='v', capsize=3, color='#252525', markersize=5)
    axl1.text(x=0.1, y=2+0.25, ha='left', va='top', size=8,
                s=('Model 1: Adjusted for difference in depth, extraction '
                   'plate,and after travel outside North America or Europe'))
    axl1.text(x=0.1, y=0.8, ha='left', va='center', size=8,
            s=('Model 2: Model 1 plus adjustment for age, sex, and vegetarian'
               ' diet'))
    axl1.text(x=0.1, y=0, ha='left', va='center', size=8, s=('p ≥ 0.05'))
    axl1.text(x=0.35, y=0, ha='left', va='center', size=8, s=('p < 0.05'))
    axl1.set_ylim(-0.5, 2.5)
    axl1.set_xlim(0, 1)
    axl1.xaxis.set_tick_params(**hide_axes)
    axl1.yaxis.set_tick_params(**hide_axes)
    
    axl2 = figl.add_subplot(gsl[1], facecolor='None', sharey=axl1)
    axl2.xaxis.set_tick_params(**hide_axes)
    axl2.yaxis.set_tick_params(**hide_axes)
    axl2.fill_between(x=[0, 1], y1=0.5, y2=1.5, color='w', edgecolor='None')
    axl2.fill_between(x=[0, 1], y1=-0.5, y2=0.5, color='#f1f1f1', edgecolor='None')
    axl2.text(0.5, 1, '7-30 days', ha='center', va='center', size=8)
    axl2.text(0.5, 0, '≥90 days', ha='center', va='center', size=8)
    axl2.text(0.5, 2, 'Time period', ha='center', va='center', size=9, weight='bold')
    axl2.set_xlim(0, 1)


def plot_alpha_curve_row(figm1, metric, sub_data, label='A', line_kws=None):
    """
    Plots the row of the alpha diveristy showing the recovery curve
    """
    # Summary parameters
    data_min = sub_data.loc[sub_data['timing'] >= 5].drop(columns=['timing'])
    data_min = np.floor(data_min.min().min() / 5) * 5
    data_max = sub_data.loc[sub_data['timing'] >= 7].drop(columns=['timing'])
    data_max = np.floor(data_max.max().max() / 5) * 5
    
    # Sets up the figure grids
    gs_m1 = figm1.add_gridspec(2, 6, height_ratios=[6, 1])
    figm11 = figm1.subfigures(1,1, facecolor='None')
    gsm11 = figm11.add_gridspec(5, 13)
    
    ylabel = \
        '\n'.join([clean_metrics.get(metric, metric.title()), '(% baseline)'])
             

    # Sets up and plots the axes
    ax_ins = []

    for abx in range(6):
        if abx == 0:
            ax0 = figm1.add_subplot(gs_m1[0, abx])
            ax1 = figm11.add_subplot(gsm11[2:4, 2 * abx + 2])
            ax0.set_ylabel(ylabel, size=10)
            ax = ax0
            axs = ax1
        else:
            ax = figm1.add_subplot(gs_m1[0, abx], sharex=ax0, sharey=ax0)
            axs = figm11.add_subplot(gsm11[2:4, 2*abx + 2], sharex=ax1, 
                                     sharey=ax1)
            _plot_trace_ax(ax, '0', sub_data, line_kws=line_kws)
            _plot_trace_ax(axs, '0', sub_data, line_kws=line_kws)
        ax.set_title(abx_order[abx], size=8)
        # antibiotic
        _plot_trace_ax(ax, str(abx), sub_data, line_kws=line_kws)
        _plot_trace_ax(axs, str(abx), sub_data, line_kws=line_kws)

        ax.set_ylim(data_min, data_max)
        axs.set_yticks(ax.get_yticks())
        axs.set_ylim(data_min, data_max)
        ax.set_xlim(-15, 365 + 15)
        ax.set_xticks(np.arange(0, 365, 30), minor=True)
        ax.set_xticks(np.arange(0, 365, 180))
        axs.set_xlim(-5, 50)
        axs.set_xticks(np.arange(0, 365, 30), minor=True)
        axs.set_xlim(-5, 50)
        
        ax.fill_between([-5, 50], data_min, data_max, facecolor='None', 
                        edgecolor='k', linestyle='-')
        ax.yaxis.set_tick_params(left=True, labelleft=(abx in {0, '0'}),
                                 right=False, labelright=False, 
                                 labelsize=8)
        axs.yaxis.set_tick_params(left=True, labelleft=False, labelsize=0)
        axs.xaxis.set_tick_params(bottom=False, labelbottom=False, 
                                  labelsize=0, which='major')
        axs.xaxis.set_tick_params(bottom=True, labelbottom=False, 
                                  labelsize=0, which='minor')
        ax.xaxis.set_tick_params(bottom=True, labelbottom=True, 
                                  labelsize=8, which='major')
        ax.xaxis.set_tick_params(bottom=True, labelbottom=False, 
                                  labelsize=0, which='minor')

        ax_ins.append(axs)
        

    inset_positions = [(0.140,  0.31, 0.075, 0.275), 
                       (0.291,  0.31, 0.075, 0.275),
                       (0.443,  0.31, 0.075, 0.275),
                       (0.595,  0.31, 0.075, 0.275),
                       (0.7455, 0.31, 0.075, 0.275),
                       (0.8955, 0.31, 0.075, 0.275),
                       ]
    for ax, pos in zip(*(ax_ins, inset_positions)):
        ax.set_position(pos)

    figm111 = figm11.subfigures(1,1, facecolor='None')
    ax0 = figm111.add_subplot(1,1,1, facecolor='None')
    ax0.xaxis.set_tick_params(bottom=False, labelbottom=False, top=False, 
                              labeltop=False) 
    ax0.yaxis.set_tick_params(bottom=False, labelbottom=False, top=False, 
                              labeltop=False) 
    ax0.set_xlim(0, 1)
    ax0.set_ylim(0, 1)
    ax0.text(0.53, -0.05, 'Days since start of treatment',
             ha='center', va='bottom', size=10)
    ax0.text(0.04, 1, label, size=12, va='top',)
    sn.despine(ax=ax0, left=True, right=True, top=True, bottom=True)

