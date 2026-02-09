import numpy as np
import pandas as pd

from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sn

# Sets up the plots so they're in a sans serif font and save correctly
# in illustrator
rcParams['font.sans-serif'] = ['Helvetica', 'Arial']
rcParams['pdf.fonttype'] = 42

# Sets up constants
hide_axes = dict(left=False, right=False, 
                 labelleft=False, labelright=False, 
                 length=0, labelsize=0)

abx_name = {'0': 'untreated',
            '1': 'nitrofurantoin', 
            '2': 'pivmecillinam',
            '3': 'co-trimoxazole',
            '4': 'ceftibuten',
            '5': 'ciprofloxacin'}

abx_color = {'0': '#969696', 
             '1': '#009E73',
             '2': '#9D654C', 
             '3': '#4E7705', 
             '4': '#098BD9',
             '5': '#CC79A7'}

def plot_sample_consort(fig, keep_count, drop_count, marker_ref=None, 
                        wrap_kws=None, drop_kws=None):
    """
    Plots a consort diagram panel for sample loss

    Parameters
    ----------
    fig: plt.Figure
        The figure where the data should be plotted
    keep_count: pd.DataFrame
        The number of samples which are dropped at each step
    drop_count: pd.DataFrame
    marker_ref: dict, pd.Series
    wrap_kws: dict
    drop_kws: dict
    """
    # Cleans up athe axes
    ax = fig.add_subplot(1,1,1, facecolor='None')
    ax.set_position((0, 0, 1, 1))
    ax.xaxis.set_tick_params(**hide_axes)
    ax.yaxis.set_tick_params(**hide_axes)
    ax.set_xlim(0, 1)
    ax.set_ylim(0.25, 8.25)
    sn.despine(ax=ax, left=True, right=True, top=True, bottom=True)

    # Pulls out the subset of the data 
    keep_small = keep_count[['step_name', 'tidy']].reset_index()
    drop_small = drop_count[['label', 'tidy', 'samples_bl']].reset_index()

    if wrap_kws is None:
        wrap_kws = dict()
    if drop_kws is None:
        drop_kws = dict()
    else:
        drop_kws.pop('marker_ref')

    # Plots the inclusions
    for idx, [step, step_name, tidy] in keep_small.iterrows():
         _build_step_wrap(ax, idx, step, step_name, tidy, **wrap_kws)
    # Plots the exclusions
    for idx, [step, label, tidy, bl] in drop_small.iterrows():
        _build_drop_step(ax, idx, step, label, tidy, bl, marker_ref=marker_ref, 
                         **drop_kws)


def build_drop_grid(fig, abx_code, abx_data):
    """
    Build a gridmap of sample loss by treatment group
    """

    # Sets up the axes
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 6, 1])
    axl = fig.add_subplot(gs[0])
    axm = fig.add_subplot(gs[1])
    axo = fig.add_subplot(gs[2], sharex=axm)

    # num_subj = len(abx_data.groupby('host_subject_id').first())
    num_subj = 15
    
    # Plots the main data
    for (mark, fc), df_ in abx_data.groupby(['marker', 'facecolor']):
        axm.scatter(x='x', y='y', s='s', 
                    color='color', 
                    marker=mark, 
                    facecolor=fc, 
                    data=df_,)
    _build_count_yaxis(axm)
    _build_count_xaxis(axm, num_subj)
    sn.despine(ax=axm, left=True, right=True, top=True, bottom=True)
    
    # Plots the by subject drops
    _build_outcome_by_subject(axo, abx_data)
    
    # Adds the number of samples at each timepoint
    _label_counts(axm, abx_data)
    
    # Adds a pretty header for the data 
    _plot_pretty_label(axl, abx_code, abx_data)


def build_grid_legend(fig):
    """
    Makes a legend for the scatter
    """
    legend_ref = pd.DataFrame(
        data=[[0, 2, 'X', 'k'],
              [0, 1, 'X', 'w'],
              [1, 2, 'd', 'k'],
              [1, 1, 'd', 'w'],
              [2, 2, 'v', 'k'],
              [2, 1, 'v', 'w'],
              [3, 2, '<', 'k'],
              [3, 1, '<', 'w'],
              [4, 2, '*', 'k'],
              [5, 2, '3', 'k']],
        columns=['x', 'y', 'marker', 'facecolor'])
    xlabels = ['Not collected\n', 
               'Not sequenced\n', 
               '>1 million\nreads',
               'Mislabeled\n', 
               'After second\nantibiotic', 
               'Extraction split\n']
    ylabels = ['Dropped', 'Participant dropped']
    
    # Sets up the legend axes
    gs = fig.add_gridspec(2, 1, height_ratios=[2, 1])
    axl1 = fig.add_subplot(gs[0])
    axl2 = fig.add_subplot(gs[1], sharex=axl1)

    # Sets up the main legend axis
    axl1.set_xlim(-0.5, 5.5)
    axl1.set_xticks(np.arange(6))
    axl1.xaxis.set_tick_params(top=False, bottom=False, labeltop=True,
                               labelbottom=False, labelsize=8, length=0)
    axl1.set_xticklabels(xlabels)

    axl1.set_ylim(0.5, 2.5)
    axl1.set_yticks([2, 1])
    axl1.yaxis.set_tick_params(**hide_axes)
    axl1.yaxis.set_tick_params(labelleft=True, labelsize=8, left=False,
                               right=False, labelright=False, length=0)
    axl1.set_yticklabels(ylabels)

    for marker, df in legend_ref.dropna().groupby('marker'):
        axl1.scatter(x='x', y='y', color='k', facecolor='facecolor',
                   data=df,marker=marker, s=48)
    
    axl2.set_ylim(0, 1)
    axl2.xaxis.set_tick_params(**hide_axes)
    axl2.yaxis.set_tick_params(**hide_axes)
    axl2.scatter(x=0, y=0.5, marker='.', color='k', s=10)
    axl2.text(0.125, 0.5, 'Sample survived', size=8, ha='left', 
              va='center')


def _build_step_wrap(ax, idx, step, step_name, tidy, 
                     stop=False, facecolor='w', 
                     edgecolor='k', top_size=10, bottom_size=8, top_color='k', 
                     bottom_color='k', linewidth=3, arrowsize=60):
    """
    A convenience function to plot the cohort drop step
    """
    y = 7.8 - idx
    ax.fill_between(x=[0.05, 0.95], 
                    y1 = y + 0.2, 
                    y2 = y - 0.2, 
                    facecolor=facecolor, 
                    edgecolor=edgecolor, 
                    linewidth=linewidth,
                    )
    ax.text(x=0.5, y=y + 0.15, s=step_name, 
            ha='center', 
            va='top', 
            size=top_size)
    ax.text(x=0.5, y=y - 0.15, s=tidy, 
            ha='center', 
            va='bottom', 
            size=bottom_size)

    if stop == False:
        ax.plot([0.1, 0.1], 
                [y - 0.2, y - 0.8], 
                color=edgecolor, 
                linewidth=linewidth)
        ax.scatter(x=0.1, 
                   y=y - 0.77, 
                   marker='v', 
                   c=edgecolor, 
                   s=arrowsize)


def _build_drop_step(ax, idx, step, label, tidy, bl, edgecolor='k', 
                     labelsize=8, 
                     marker_ref=None, drop_facecolor='#dbdbdb', 
                     linewidth=3,
                     arrowsize=60):
    """
    Plots the drop out box
    """
    y = 7.25 - idx

    # Plots the connection and call out
    ax.plot([0.1, 0.15], [y, y], 
            color=edgecolor, 
            linewidth=linewidth)
    ax.scatter(x=[0.145], 
               y=[y], 
               marker='>', 
               color=edgecolor, 
               s=arrowsize, 
               zorder=3,)
    ax.fill_between(x=[0.165, 0.95], 
                    y1= y + 0.125, 
                    y2=y - 0.125, 
                    facecolor=drop_facecolor)

    # Sets up the drop out marker
    if marker_ref is None:
        marker = 'None'
    else:
        marker = marker_ref.get(step, 'None')

    # Plots the drop out scatter
    ax.scatter(x=0.2,
               y=y + 0.01, 
               marker=marker, 
               color='k')
    ax.text(x=0.225, y=y, s=tidy, 
            ha='left', va='center', size=labelsize)
    if bl > 0:
        ax.scatter(x=0.66, 
                   y=y + 0.01, 
                   marker = marker,
                   color='k', 
                   facecolor='w')
        ax.text(x=0.69, y=y, s=f'{bl} samples', 
                ha='left', va='center', size=labelsize)


def _build_count_yaxis(ax):
    """
    Cleans up the yaxis for the antibiotic heatmap
    """
    # Sets up the time frame ticks
    ax.yaxis.set_tick_params(**hide_axes)
    ax.yaxis.set_tick_params(labelleft=True, labelsize=8)
    ax.set_yticks(np.arange(7), minor=False)
    ax.set_yticklabels(['-1', '+7', '+14', '+30', "+90", '+181', '+365'])

    # Sets up the minor ticks to create the grid
    ax.yaxis.set_tick_params(**hide_axes, gridOn=True, which='minor')
    ax.set_yticks(np.arange(8) - 0.5, minor=True)

    ax.set_ylim(6.6, -0.6)

    ax.set_ylabel('Timepoint', size=10)


def _build_count_xaxis(ax, num_subj):
    """
    Cleans up the xaxis for the anitbiotic heatmap
    """
    # Hides the major ticks
    ax.xaxis.set_tick_params(**hide_axes)

    # Sets up the major ticks for the grid
    ax.set_xticks(np.arange(num_subj + 3) - 0.5, minor=True)
    ax.xaxis.set_tick_params(**hide_axes, gridOn=True, which='minor')
    ax.set_xlim(-0.6, num_subj - 0.3)

    
    # Adds line indicating the baseline sample
    ax.plot([-1, num_subj + 1], [0.5, 0.5], 'k-', linewidth=2)


def _build_outcome_by_subject(ax, abx_data):
    """
    Counts and labels participant drops
    """
    # Determines how many samples each subject has, or what step they 
    # were dropped at
    drop_count = \
        abx_data.groupby(['host_subject_id', 'drop_count'], sort=False)
    drop_count = drop_count[['marker', 'x']].first()
    drop_count.reset_index(level='drop_count', inplace=True)
    drop_count['drop_count'] = 7 - drop_count['drop_count'].astype(int)
    drop_count['mark'] = drop_count['drop_count'] == 0

    for mark, plot_df in drop_count.groupby('marker'):
        if mark == '.':
            for _, [x, c] in plot_df[['x', 'drop_count']].iterrows():
                ax.text(x=x, y=0.5, s=c, ha='center', va='center', size=8)
        else:
            ax.scatter(x=plot_df['x'], y=plot_df['x'] * 0 + 0.5, 
                       marker=mark, color='k', s=36)

    ax.xaxis.set_tick_params(**hide_axes)
    ax.yaxis.set_tick_params(**hide_axes)
    ax.set_yticks([-0.1, 1.1], minor=True)
    ax.xaxis.set_tick_params(**hide_axes, gridOn=False, which='minor')
    ax.yaxis.set_tick_params(**hide_axes, gridOn=False, which='minor')

    sn.despine(ax=ax, left=True, right=True, top=True, bottom=True)


def _label_counts(ax, abx_data):
    """
    Adds a hidden axis identifying the number of samples in the dataset
    """
    # Determines hte number of samples at each timepoint
    time_count = abx_data.set_index('timepoint')['marker'].isin(['.'])
    time_count = time_count.groupby('timepoint').sum()

    axt = ax.twinx()
    axt.set_ylim(ax.get_ylim())
    axt.set_yticks(ax.get_yticks())
    axt.yaxis.set_tick_params(**hide_axes)
    axt.yaxis.set_tick_params(labelright=True, labelsize=8)
    axt.set_yticklabels(time_count.values)

    sn.despine(ax=axt, left=True, right=True, top=True, bottom=True)

def _plot_pretty_label(ax, abx, abx_data):
    """
    Makes a pretty title plot of the data
    """

    final_count = abx_data.groupby('host_subject_id')['marker'].first()
    final_count = final_count.isin(['.']).sum()
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.xaxis.set_tick_params(**hide_axes)
    ax.yaxis.set_tick_params(**hide_axes)
    sn.despine(ax=ax, left=True, right=True, bottom=True, top=True)

    ax.fill_between([0, 1], 0, 1, color=abx_color[abx], alpha=1/3)
    ax.text(x=0.5, y=0.5, 
            s='{0} (n={1})'.format(abx_name[abx].title(), final_count),
            ha='center', va='center', size=10)

       