#!/usr/bin/env python
# coding=utf-8
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import re
import argparse
import matplotlib as mpl
mpl.use('agg')

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--plot_type', type=str, default=64, help='DimPlot or FeaturePlot')
parser.add_argument('-i', '--plot_item', type=str, default=64, help='plot item')
parser.add_argument('--plot_order', type=str, default=64, help='plot order')
parser.add_argument('-o', '--output', type=str, default=64, help='output prefix')

args = parser.parse_args()
plot_type = args.plot_type
plot_item = args.plot_item
plot_order = args.plot_order
prefix = args.output

if plot_type == 'DimPlot':
    plot_order = plot_order.split(' ')
    print(plot_order)
mycolor = ["#B2D4F9","#AF896E","#6AB070","#eaeaea"] #Replace as needed


cm = plt.cm.Spectral_r

df_meta = pd.read_csv('df_meta_tmp.txt', index_col=0, sep='\t')
df_meta = df_meta.drop(['row', 'col'], axis=1)
df_pos = pd.read_csv('df_pos_tmp.txt', index_col=0, sep='\t')

df_pos = df_pos[['row', 'col']]
df_pos['row'] = df_pos['row'] - min(df_pos['row']) + 1          
df_pos['col'] = df_pos['col'] - min(df_pos['col']) + 1

# Swap rows and columns 
df_pos['tmp'] = df_pos['row'].copy()
df_pos['row'] = df_pos['col']
df_pos['col'] = df_pos['tmp']
df_pos = df_pos.drop('tmp', axis=1)    
df_pos['col'] = max(df_pos['col']) + 1 - df_pos['col']


spot_row = max(df_pos['row']) - min(df_pos['row']) + 1        
spot_col = max(df_pos['col']) - min(df_pos['col']) + 1
max_spot = max(spot_row, spot_col)

df_plot = pd.merge(df_pos, df_meta, left_index=True, right_index=True)

Sizef = 5
# featureplot
def featureplot_inhouse(df_plot, item, prefix):
    fig, ax = plt.subplots(figsize=(Sizef * spot_row / spot_col, Sizef))
    ax.set_xlim((0 ,spot_row+1))
    ax.set_ylim((0, spot_col+1))
    r = 0.55
    r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
    marker_size = np.pi * r_**2
    ax.scatter(x=df_plot['row'], y=df_plot['col'], c=list(df_plot[item]), s=marker_size, edgecolors='black', linewidth=0, cmap=cm)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.gca().add_patch(plt.Rectangle(xy=(0.5,0.5),
                                  width=spot_row,
                                  height=spot_col,
                                  edgecolor='#bdbdbd',
                                  fill=False, linewidth=1.25))
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0.1, wspace=0.1)
    ax.figure.savefig(f'{prefix}.png',dpi=300)


def dimplot_inhouse(df_plot, item, prefix, plot_order, cols=mycolor):
    cur_dic = {plot_order[idx]: cols[idx] for idx in range(len(plot_order))}
    print(cur_dic)
    fig, ax = plt.subplots(figsize=(Sizef * spot_row / spot_col, Sizef))
    ax.set_xlim((0 ,spot_row+1))
    ax.set_ylim((0, spot_col+1))
    r = 0.55
    r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
    marker_size = np.pi * r_**2

    # Plot the desired cluster
    for cur_item in cur_dic.keys():
        df_highlight = df_plot[df_plot[item] == cur_item]
        ax.scatter(x=df_highlight['row'], y=df_highlight['col'], c=cur_dic[cur_item], s=marker_size, edgecolors='none')
    
    # Plot all other points in gray
    back_ind = df_plot.copy()
    for cur_item in cur_dic.keys():
        back_ind = back_ind[back_ind[item] != cur_item]
    ax.scatter(x=back_ind['row'], y=back_ind['col'], c="#eaeaea", s=marker_size, edgecolors='none')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.gca().add_patch(plt.Rectangle(xy=(0.5,0.5),
                                    width=spot_row,
                                    height=spot_col,
                                    edgecolor='#bdbdbd',
                                    fill=False, linewidth=1.25))
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0.1, wspace=0.1)
    ax.figure.savefig(f'{prefix}.png',dpi=200 )

    fig, ax = plt.subplots(figsize=(Sizef * spot_row / spot_col, Sizef))
    ax.set_xlim((0 ,spot_row+1))
    ax.set_ylim((0, spot_col+1))
    r = 0.55
    r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
    marker_size = np.pi * r_**2

    for cur_item in cur_dic.keys():
        df_highlight = df_plot[df_plot[item] == cur_item]
        ax.scatter(x=df_highlight['row'], y=df_highlight['col'], c=cur_dic[cur_item], s=marker_size, edgecolors='none', label=cur_item)
        # plot all other points in gray
    back_ind = df_plot.copy()
    for cur_item in cur_dic.keys():
        back_ind = back_ind[back_ind[item] != cur_item]
    ax.scatter(x=back_ind['row'], y=back_ind['col'], c="#a1a499", s=marker_size, edgecolors='none')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.gca().add_patch(plt.Rectangle(xy=(0.5,0.5),
                                    width=spot_row,
                                    height=spot_col,
                                    edgecolor='#bdbdbd',
                                    fill=False, linewidth=1.25))
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0.1, wspace=0.1)
    plt.legend(prop={'family':'SimHei','size':6},markerscale=4)
    ax.figure.savefig(f'{prefix}.legend.png',dpi=200)


if plot_type == 'DimPlot':
    df_plot[plot_item] = df_plot[plot_item].astype('str')
    dimplot_inhouse(df_plot, plot_item, prefix, plot_order)
elif plot_type == 'FeaturePlot':
    featureplot_inhouse(df_plot, plot_item, prefix)
