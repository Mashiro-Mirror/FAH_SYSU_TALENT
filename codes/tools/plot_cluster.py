#!/usr/bin/env python
# coding=utf-8
import os
import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import re
import argparse
mpl.use('agg')

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--plot_type', type=str, default=64, help='DimPlot or FeaturePlot')
parser.add_argument('-i', '--plot_item', type=str, default=64, help='plot item')
parser.add_argument('--plot_order', type=str, default=64, help='plot order')
parser.add_argument('-o', '--output', type=str, default=64, help='output prefix')
parser.add_argument('-n', '--name', type=str, default=64, help='output name')
parser.add_argument('--plot_color', type=str, default=64, help='color')
parser.add_argument('--plot_tls_contour', type=str, default=64, help='color')

args = parser.parse_args()
plot_type = args.plot_type
plot_item = args.plot_item
plot_order = args.plot_order
mycolor = args.plot_color
plot_tls_contour = args.plot_tls_contour

name = args.name

prefix = args.output

namels = prefix.split("/")
sampid = namels[-1].split("_")[0]

if plot_type == 'DimPlot':
    plot_order = plot_order.split(' ')
    print(plot_order)
    mycolor = mycolor.split(' ')
    print(mycolor)
    
col16 = ["#EE3743", "#33a02c", "#0B7FAB", "#6a3d9a", "#ff5084", "#ff7f00", "#01ceff", "#F6A9BD", "#fdbf6f", "#c51b7d",  
                            "#6F6C9E", "#01665e", "#a6cee3", "#b2df8a", "#ffff99", "#bf812d",  "#999999", "#BB7DB2", "#824615", "#ffff33", "#86DBD4",
                            "#BFE2E3","#A1CFFA","#78BDAD","#D45651","#397A7F","#F0918E","#EEE8DA","#1F5392","#A0BFAF",
                            "#AE98D6","#ECCBDC","#54BAD3","#8b4a4b","#DB896C","#AABAC2","#ffae3b",
                            '#CCCCCC','#B5B5B5','#A092C3','#DDA0DD','#03A4C6','#7AC5CD','#A3D9ED','#00E5EE','#F9EBDA','#F5DEB3',
                            '#98689E','#E84115','#FFB5C5','#00FF7F','#F9D01C','#B03060','#00ABDC','#D2691E','#03A464','#FF7F00',
                            '#8968CD','#1C5B75']

# cm = plt.cm.Spectral_r
cm = mpl.colors.LinearSegmentedColormap.from_list('custom_grey_red', ['#d3d3d3', '#ff0000'])

df_meta = pd.read_csv(f'{prefix}/df_meta_tmp.txt', index_col=0, sep='\t')
df_meta = df_meta.drop(['row', 'col'], axis=1)
df_pos = pd.read_csv(f'{prefix}/df_pos_tmp.txt', index_col=0, sep='\t')

df_pos = df_pos[['row', 'col']]
df_pos['row'] = df_pos['row'] - min(df_pos['row']) + 1          # 预留5格作为边框
df_pos['col'] = df_pos['col'] - min(df_pos['col']) + 1

# 将行列交换 
df_pos['tmp'] = df_pos['row'].copy()
df_pos['row'] = df_pos['col']
df_pos['col'] = df_pos['tmp']
df_pos = df_pos.drop('tmp', axis=1)     # 删除tmp列
df_pos['col'] = max(df_pos['col']) + 1 - df_pos['col']


spot_row = max(df_pos['row']) - min(df_pos['row']) + 1        # 预留5格作为边框
spot_col = max(df_pos['col']) - min(df_pos['col']) + 1
max_spot = max(spot_row, spot_col)
df_plot = pd.merge(df_pos, df_meta, left_index=True, right_index=True)


bin_f = np.full((spot_col+2, spot_row+2), 0, dtype="uint8")
bin_f[df_pos['col'].values, df_pos['row'].values] = 1
kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))       # 构建形态学算子
erode = cv2.dilate(bin_f, kernel)    # 腐蚀图像
edge = erode - bin_f                # 获取边界
edgeXY = np.where(edge == 1)      # 获取edge的坐标


Sizef = 5
backcolor = "#eaeaea"
# featureplot

def plot_contour(df_plot, spot_row, spot_col, drawtype, lineC='#000000',lineW=1, ax = None):
    """
        df_plot:输入完整的包含待勾画轮廓部分的数据
        spot_row, spot_col:与输入数据对应的spot行和列
        drawtype:勾画轮廓的类型,如 slide:整张slide/TLS:所有TLS或单个TLS
        color='#000000',linewidth=1 绘图参数:轮廓线颜色与宽度
    """ 
    def drawedge(df_data, spot_row, spot_col, lineC='#000000',lineW=1, ax = None):
        # 构建灰度图像
        bin_f = np.full((spot_col+2, spot_row+2), 0, dtype="uint8")         # 构建图像  
        bin_f[df_data['col'].values, df_data['row'].values] = 1             
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))       # 构建形态学算子
        erode = cv2.dilate(bin_f, kernel)                                   # 扩增图像（为了勾画TLS外一圈的轮廓）

        # 识别灰度图像轮廓坐标
        _, thresh = cv2.threshold(erode,  0.8, 1, cv2.THRESH_BINARY)
        plot_contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)      # 4. 执行边缘检测，获取轮廓
        
        for contour in plot_contours:
            # 取出连续相邻的坐标进行绘制
            curlist = [[i, i+1 if i+1 < len(contour) else 0] for i in range(len(contour))]
            plot_x = [[contour[cur_2][0][0]+1, contour[cur_1][0][0]+1] for cur_1, cur_2 in curlist]
            plot_y = [[contour[cur_2][0][1]+1, contour[cur_1][0][1]+1] for cur_1, cur_2 in curlist]
            ax.plot(plot_x, plot_y, color=lineC,linewidth=lineW)

    if drawtype == 'TLS':
        df_data = df_plot[~df_plot['TLS'].str.endswith('_NA')]          # 取出对应数据
        tlslist = df_data.TLS.unique()
        for ti in tlslist:
            subda = df_plot.loc[df_plot.TLS == ti]
            drawedge(subda, spot_row, spot_col, lineC, lineW, ax = ax)
    else:
        df_data = df_plot
        drawedge(df_data, spot_row, spot_col, lineC, lineW, ax = ax)


def featureplot_inhouse(df_plot, item, prefix):
    vmin = df_plot[item].min()
    vmax = df_plot[item].max()
    vcenter = np.percentile(df_plot[item], 25)
    if not (vmin < vcenter < vmax):
        vcenter = (vmin + vmax) / 2
    norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)  # Adjust vcenter to lower quartile for better color contrast
      # Adjust vcenter to lower quartile for better color contrast  # Adjust color scale to make lower values lighter and higher values more vibrant
    # fig, ax = plt.subplots(figsize=(Sizef * spot_row / spot_col, Sizef))
    # ax.set_xlim((0 ,spot_row+3))
    # ax.set_ylim((0, spot_col+3))
    # r = 0.55
    # r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
    # marker_size = np.pi * r_**2
    # ax.scatter(x=df_plot['row']+1, y=df_plot['col']+1, c=list(df_plot[item]), s=marker_size, edgecolors='none',                         linewidth=0, cmap=cm)
    # ax.scatter(x=edgeXY[1]+1, y=edgeXY[0]+1, c="#000000", s=marker_size, edgecolors="none")
    # ax.set_xticks([])
    # ax.set_yticks([])
    # ax.axis('off')
    # plt.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0.1, wspace=0.1)
    # # ax.figure.savefig(f'{prefix}.png', dpi=800)
    # ax.figure.savefig(f'{prefix}/{name}.pdf')
    
    # 主要得到spot_row的信息
    
    fig, ax = plt.subplots(figsize=(Sizef * spot_row / spot_col, Sizef))
    ax.set_xlim((0 ,spot_row+3))
    ax.set_ylim((0, spot_col+3))
    r = 0.55
    r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
    marker_size = np.pi * r_**2
    ax.scatter(x=df_plot['row']+1, y=df_plot['col']+1, c=list(df_plot[item]), s=marker_size, edgecolors='none', linewidth=0, cmap=cm, alpha=0.8, norm=norm)
    if plot_tls_contour == "T":
        plot_contour(df_plot, spot_row, spot_col, 'TLS', '#000000', 0.5, ax = ax)
    # ax.scatter(x=edgeXY[1]+1, y=edgeXY[0]+1, c="#000000", s=marker_size, edgecolors="none")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0.1, wspace=0.1)
    ax.figure.savefig(f'{prefix}/{name}.png', dpi=400)
    ax.figure.savefig(f'{prefix}/{name}.pdf')




def dimplot_inhouse(df_plot, item, prefix, plot_order, cols=mycolor):
    cur_dic = {plot_order[idx]: cols[int(idx)] for idx in range(len(plot_order))}
    print(cur_dic)
    fig, ax = plt.subplots(figsize=(Sizef * spot_row / spot_col, Sizef))
    ax.set_xlim((0 ,spot_row+3))
    ax.set_ylim((0, spot_col+3))
    r = 0.55
    r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
    marker_size = np.pi * r_**2

    # 首先将想要画的cluster绘制上去
    for cur_item in cur_dic.keys():
        df_highlight = df_plot[df_plot[item] == cur_item]
        ax.scatter(x=df_highlight['row']+1, y=df_highlight['col']+1, c=cur_dic[cur_item], s=marker_size, edgecolors='none')
    
    # 然后将其他的点全部绘制成为灰色
    back_ind = df_plot.copy()
    for cur_item in cur_dic.keys():
        back_ind = back_ind[back_ind[item] != cur_item]
    ax.scatter(x=back_ind['row']+1, y=back_ind['col']+1, c=backcolor, s=marker_size, edgecolors='none')
    if plot_tls_contour == "T":
        plot_contour(df_plot, spot_row, spot_col, 'TLS', '#000000', 0.5, ax = ax)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0.1, wspace=0.1)
    # ax.figure.savefig(f'{prefix}/{name}.png', dpi=1200)
    ax.figure.savefig(f'{prefix}/{name}.pdf')
    # ax.figure.savefig(f'{prefix}/{name}.png', dpi=800)

    # ax.figure.savefig(f'{prefix}/{name}.png', dpi=800)

    fig, ax = plt.subplots(figsize=(Sizef * spot_row / spot_col, Sizef))
    ax.set_xlim((0 ,spot_row+3))
    ax.set_ylim((0, spot_col+3))
    r = 0.55
    r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
    marker_size = np.pi * r_**2

    for cur_item in cur_dic.keys():
        df_highlight = df_plot[df_plot[item] == cur_item]
        ax.scatter(x=df_highlight['row']+1, y=df_highlight['col']+1, c=cur_dic[cur_item], s=marker_size, edgecolors='none', label=cur_item)
        # 然后将其他的点全部绘制成为灰色
    back_ind = df_plot.copy()
    for cur_item in cur_dic.keys():
        back_ind = back_ind[back_ind[item] != cur_item]
    ax.scatter(x=back_ind['row']+1, y=back_ind['col']+1, c=backcolor, s=marker_size, edgecolors='none')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0.1, wspace=0.1)
    plt.legend(prop={'family':'sans','size':6},markerscale=4)
    ax.figure.savefig(f'{prefix}/{name}.legend.png', dpi=400)
    ax.figure.savefig(f'{prefix}/{name}.legend.pdf')


def plot_colorbar(df_plot, item, prefix):
    fig, ax = plt.subplots(figsize=(Sizef * spot_row / spot_col, Sizef))
    ax.set_xlim((0, spot_row+1))
    ax.set_ylim((0, spot_col+1))
    r = 0.55
    r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
    marker_size = np.pi * r_**2
    sc = ax.scatter(x=df_plot['row'], y=df_plot['col'], c=list(df_plot[item]), s=marker_size, edgecolors='black', linewidth=0, cmap=cm)
    plt.cla()                                                           # 将绘制的图都清除

    # formatter = mpl.ticker.StrMethodFormatter('{x:.0f}')                # 设置为刻度为整数
    # cbar = fig.colorbar(sc, shrink=0.8, orientation='horizontal', ax=ax, format = formatter)
    # cbar.set_ticks(np.arange(int(cbar.vmin), int(cbar.vmax+1)+1, (int(cbar.vmax+1) - int(cbar.vmin))/4))      # 设置整数刻度
    cbar = fig.colorbar(sc, shrink=0.8, orientation='horizontal', ax=ax)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    plt.subplots_adjust(left=0, bottom=0.4, right=1, top=1, hspace=0.1, wspace=0.1)
    plt.legend(prop={'family': 'sans', 'size': 6}, markerscale=4)
    ax.figure.savefig(f'{prefix}/{name}.colorbar.png', dpi=800)
    ax.figure.savefig(f'{prefix}/{name}.colorbar.pdf')



if plot_type == 'DimPlot':
    df_plot[plot_item] = df_plot[plot_item].astype('str')
    dimplot_inhouse(df_plot, plot_item, prefix, plot_order)
elif plot_type == 'FeaturePlot':
    featureplot_inhouse(df_plot, plot_item, prefix)
elif plot_type == 'Colorbar':
    plot_colorbar(df_plot, plot_item, prefix)
