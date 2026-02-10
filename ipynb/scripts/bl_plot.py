import numpy as np
import pandas as pd
import skbio

from matplotlib import rcParams
from matplotlib import ticker
import matplotlib.pyplot as plt
import seaborn as sn


tick_in = dict(left=True, right=True, labelleft=False, labelright=False, 
               labelsize=0, tickdir='in',)
hide_axes = dict(left=False, right=False, labelleft=False, labelright=False, 
                 labelsize=0, length=0)

def plot_ordination(fig, pc0, pc1, hue, palette, data, perc_exp, 
                    scatter_kws=None, box_kws=None):
    """
    Makes an ordination plot
    """
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 4], width_ratios=[4, 1])
    axm = fig.add_subplot(gs[1, 0])
    ax0 = fig.add_subplot(gs[0, 0], sharex=axm, facecolor='None')
    ax1 = fig.add_subplot(gs[1, 1], sharey=axm, facecolor='None')

    scatter_kwargs = dict(edgecolor='None',
                          s=2,)
    if isinstance(scatter_kws, dict):
        scatter_kwargs.update(scatter_kws)
    box_kwargs = dict(linewidth=0.5,
                      fliersize=0.5)
    if isinstance(box_kws, dict):
        box_kwargs.update(box_kws)
    
    sn.scatterplot(ax=axm, x=pc0, y=pc1, data=data, hue=hue, palette=palette, 
                   **scatter_kwargs)
    sn.boxplot(ax=ax0, x=pc0, y=hue, palette=palette, data=data, **box_kwargs)
    sn.boxplot(ax=ax1, y=pc1, x=hue, palette=palette, data=data, **box_kwargs)
    axm.get_legend().set_visible(False)
    axm.xaxis.set_major_locator(ticker.LinearLocator(5))
    axm.xaxis.set_tick_params(**tick_in)
    axm.yaxis.set_major_locator(ticker.LinearLocator(5))
    axm.yaxis.set_tick_params(**tick_in)
    axm.set_xlabel(f'PC 1 ({{{pc0}:1.1%}})'.format(**perc_exp), size=7)
    axm.set_ylabel(f'PC 2 ({{{pc1}:1.1%}})'.format(**perc_exp), size=7)
    
    for ax_ in [ax0.xaxis, ax0.yaxis, ax1.xaxis, ax1.yaxis]:
        ax_.set_tick_params(**hide_axes)
        ax_.get_label().set_visible(False)

    sn.despine(ax=ax0, left=True, right=True, top=True, bottom=True)
    sn.despine(ax=ax1, left=True, right=True, top=True, bottom=True)
