#!/usr/bin/env python
# coding=utf-8
import cv2
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import re
mpl.use('agg')
ax = plt.gca()

########## Plotting whole slide of ST data with highlight of TLS region ##########
def plot_contour(df_plot, spot_row, spot_col, drawtype, lineC='#000000', lineW=1):
    """
        df_plot: Input the complete dataset that includes the section to be outlined
        spot_row, spot_col: The rows and columns corresponding to the input data
        drawtype: The type of contour to be outlined, such as slide: entire slide/TLS: all TLS or individual TLS
        color='#000000',linewidth=1 Plotting parameters: contour line color and width
    """

    def drawedge(df_data, spot_row, spot_col, lineC='#000000', lineW=1):
        # Constructing a grayscale image
        bin_f = np.full((spot_col + 2, spot_row + 2), 0, dtype="uint8")  
        bin_f[df_data['col'].values, df_data['row'].values] = 1
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))  
        erode = cv2.dilate(bin_f, kernel)  

        # Identifying coordinates of grayscale image contours
        _, thresh = cv2.threshold(erode, 0.8, 1, cv2.THRESH_BINARY)
        plot_contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)  

        for contour in plot_contours:
            # Plotting continuous adjacent coordinates
            curlist = [[i, i + 1 if i + 1 < len(contour) else 0] for i in range(len(contour))]
            plot_x = [[contour[cur_2][0][0] + 1, contour[cur_1][0][0] + 1] for cur_1, cur_2 in curlist]
            plot_y = [[contour[cur_2][0][1] + 1, contour[cur_1][0][1] + 1] for cur_1, cur_2 in curlist]
            ax.plot(plot_x, plot_y, color=lineC, linewidth=lineW)

    if drawtype == 'TLS':
        df_data = df_plot[~df_plot['TLS'].str.endswith('_NA')]  
        tlslist = df_data.TLS.unique()
        for ti in tlslist:
            subda = df_plot.loc[df_plot.TLS == ti]
            drawedge(subda, spot_row, spot_col, lineC, lineW)
    else:
        df_data = df_plot
        drawedge(df_data, spot_row, spot_col, lineC, lineW)

def featureplot_inhouse1(df_plot, spot_row, spot_col, prefix):
    Sizef = 5
    cm = plt.cm.Spectral_r
    # Obtaining the coordinates of the boundaries
    bin_f = np.full((spot_col + 2, spot_row + 2), 0, dtype="uint8")
    bin_f[df_plot['col'].values, df_plot['row'].values] = 1
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))  
    erode = cv2.erode(bin_f, kernel)
    edge = bin_f - erode
    edgeXY = np.where(edge == 1) 

    # Obtaining the coordinates of the boundaries belong to "B"
    df_plotB = df_plot.loc[df_plot.Bin_Region == 'B']
    bin_f = np.full((spot_col + 2, spot_row + 2), 0, dtype="uint8")
    bin_f[df_plotB['col'].values, df_plotB['row'].values] = 1
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))  
    erode1 = cv2.erode(bin_f, kernel) 
    edge1 = bin_f - erode1  
    edgeXY1 = np.where(edge1 == 1) 

    sampid = prefix.split('/')[-1]
    df_plot = df_plot.loc[df_plot.TLS != sampid + "_NA"]

    # Plotting
    fig, ax = plt.subplots(figsize=(Sizef * spot_row / spot_col, Sizef))
    ax.set_xlim((0, spot_row + 3))
    ax.set_ylim((0, spot_col + 3))
    r = 0.55
    r_ = ax.transData.transform([r, 0])[0] - ax.transData.transform([0, 0])[0]
    marker_size = np.pi * r_ ** 2
    for mm in list(Mycolor.keys()):
        df_TLSda_ma = df_plot.loc[df_plot.TLS_use == mm]
        df_TLSda_ma = df_TLSda_ma.loc[df_TLSda_ma.TLS != sampid + "_NA"]
        ax.scatter(x=df_TLSda_ma['row'] + 1, y=df_TLSda_ma['col'] + 1, c=Mycolor.get(mm), s=marker_size,
                   edgecolors="none")
    ax.scatter(x=edgeXY1[1]+1, y=edgeXY1[0]+1, c="#000000", s=marker_size, edgecolors="none")
    ax.scatter(x=edgeXY[1]+1, y=edgeXY[0]+1, c="#000000", s=marker_size, edgecolors="none")
    plot_contour(df_pos, spot_row, spot_col, 'slide', '#000000', 0.5)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    # plt.gca().add_patch(plt.Rectangle(xy=(0.5,0.5),
    #                               width=spot_row+2,
    #                               height=spot_col+2,
    #                               edgecolor='#bdbdbd',
    #                               fill=False, linewidth=1.25))
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0.1, wspace=0.1)
    ax.figure.savefig(f'{prefix}.pdf')
    plt.close()


pathans = "/Users/ctu/Desktop/ST/TLS/"
HCC_meta_path = "/Users/ctu/Desktop/ST/TLS/"
samplelist = ["pt03","pt05","pt06","pt08","pt10","pt14","pt15","pt16","pt17"]

Mycolor = {'TLS': '#4DBBD5','Tumor': '#E8E5F0','Peritumor': '#F8F5EC'} 

for ii in range(len(samplelist)):
    sampid = samplelist[ii]

    df_meta = pd.read_csv(HCC_meta_path + sampid + "_metadata.txt", sep='\t', index_col=0, low_memory=False)

    df_pos = df_meta[['row', 'col', 'TLS', 'TLS_maturity', 'Bin_Region']]
    df_pos['row'] = df_pos['row'] - min(df_pos['row']) + 1  
    df_pos['col'] = df_pos['col'] - min(df_pos['col']) + 1
    df_pos['tmp'] = df_pos['row'].copy()
    df_pos['row'] = df_pos['col']
    df_pos['col'] = df_pos['tmp']
    df_pos = df_pos.drop('tmp', axis=1) 
    df_pos['col'] = max(df_pos['col']) + 1 - df_pos['col']
    max_col = max(df_pos['col'])
    spot_row = max(df_pos['row']) - min(df_pos['row']) + 1  
    spot_col = max(df_pos['col']) - min(df_pos['col']) + 1
    max_spot = max(spot_row, spot_col)
    # df_plot = pd.merge(df_pos, df_meta, left_index=True, right_index=True)
    prefix = pathans + sampid
    featureplot_inhouse1(df_pos, spot_row, spot_col, prefix=prefix)  


########## Plotting TLS region only with celltypes annotation ##########
def plot_contour(df_plot, spot_row, spot_col, drawtype, lineC='#000000',lineW=1):
    """
        df_plot: Input the complete dataset that includes the section to be outlined
        spot_row, spot_col: The rows and columns corresponding to the input data
        drawtype: The type of contour to be outlined, such as slide: entire slide/TLS: all TLS or individual TLS
        color='#000000',linewidth=1 Plotting parameters: contour line color and width
    """ 
    def drawedge(df_data, spot_row, spot_col, lineC='#000000',lineW=1):
        # Constructing a grayscale image
        bin_f = np.full((spot_col+2, spot_row+2), 0, dtype="uint8")         
        bin_f[df_data['col'].values, df_data['row'].values] = 1             
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))       
        erode = cv2.dilate(bin_f, kernel)                                   

        # Identifying coordinates of grayscale image contours
        _, thresh = cv2.threshold(erode, 0.8, 1, cv2.THRESH_BINARY)
        plot_contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)   
        
        for contour in plot_contours:
            # Plotting continuous adjacent coordinates
            curlist = [[i, i+1 if i+1 < len(contour) else 0] for i in range(len(contour))]
            plot_x = [[contour[cur_2][0][0]+1, contour[cur_1][0][0]+1] for cur_1, cur_2 in curlist]
            plot_y = [[contour[cur_2][0][1]+1, contour[cur_1][0][1]+1] for cur_1, cur_2 in curlist]
            ax.plot(plot_x, plot_y, color=lineC,linewidth=lineW)

    if drawtype == 'TLS':
        df_data = df_plot[~df_plot['TLS'].str.endswith('_NA')]         
        tlslist = df_data.TLS.unique()
        for ti in tlslist:
            subda = df_plot.loc[df_plot.TLS == ti]
            drawedge(subda, spot_row, spot_col, lineC, lineW)
    else:
        df_data = df_plot
        drawedge(df_data, spot_row, spot_col, lineC, lineW)

def adjrowcol(df_pos, TLSxyda):
    tls_tmp_row, tls_tmp_col = max(TLSxyda['row']) - min(TLSxyda['row']) + 1, max(TLSxyda['col']) - min(TLSxyda['col']) + 1
    add_tmp_row, add_tmp_col = (int(tls_tmp_row/0.618) - tls_tmp_row)//2 + 1, (int(tls_tmp_col/0.618) - tls_tmp_col)//2 + 1
    min_row, min_col = min(TLSxyda['row']) - add_tmp_row + 1, min(TLSxyda['col']) - add_tmp_col + 1
    max_row, max_col = max(TLSxyda['row']) + add_tmp_row + 1, max(TLSxyda['col']) + add_tmp_col + 1
   
    min_row = min_row if min_row > 0 else 1
    min_col = min_col if min_col > 0 else 1
    max_row = max_row if max_row < max(df_pos['row']) else max(df_pos['row'])
    max_col = max_col if max_col < max(df_pos['col']) else max(df_pos['col'])
   
    flag = (df_pos['row'] >= min_row).values &\
           (df_pos['row'] <= max_row).values &\
           (df_pos['col'] >= min_col).values &\
           (df_pos['col'] <= max_col).values

    TLSxy1 = df_pos[flag]
    TLSxy1['row'] = TLSxy1['row'] - min(TLSxy1['row']) + 1
    TLSxy1['col'] = TLSxy1['col'] - min(TLSxy1['col']) + 1
    tls_row = max(TLSxy1['row']) - min(TLSxy1['row']) + 1
    tls_col = max(TLSxy1['col']) - min(TLSxy1['col']) + 1
    return tls_row, tls_col, TLSxy1, flag


pathans = "/Users/ctu/Desktop/ST/TLS/"
HCC_meta_path = "/Users/ctu/Desktop/ST/TLS/TLS_region_only"
samplelist = ["pt08","pt10","pt15"]

Sel_CellTypeList = ['Plasma','B_cell','Other_cells']
mycolor_Plasma_only = ['#d9644a','#6BBF98','#eaeaea']

for ii in range(len(samplelist)):
    sampid  = samplelist[ii]
    os.system(f'mkdir {sampid}')
    df_meta = pd.read_csv(HCC_meta_path+sampid+".txt" ,sep='\t', index_col=0,low_memory=False)
    
    df_pos  = df_meta[['row', 'col', 'TLS', 'CellSubType_new']]
    min_row = min(df_pos['row'])
    min_col = min(df_pos['col'])
    df_pos['row'] = df_pos['row'] - min(df_pos['row']) + 1          
    df_pos['col'] = df_pos['col'] - min(df_pos['col']) + 1
    df_pos['tmp'] = df_pos['row'].copy()
    df_pos['row'] = df_pos['col']
    df_pos['col'] = df_pos['tmp']
    df_pos = df_pos.drop('tmp', axis=1)                             
    df_pos['col'] = max(df_pos['col']) + 1 - df_pos['col']
    max_col  = max(df_pos['col'])
    spot_row = max(df_pos['row']) - min(df_pos['row']) + 1          
    spot_col = max(df_pos['col']) - min(df_pos['col']) + 1
    max_spot = max(spot_row, spot_col)
    # df_plot = pd.merge(df_pos, df_meta, left_index=True, right_index=True)
    df_pos1 = df_pos.loc[df_pos.TLS != sampid+'_NA']
    TLSlist = df_pos1.TLS.unique()
    for sub_tls in TLSlist:
        TLSxy = df_pos1[df_pos1.TLS == sub_tls]
        spot_row, spot_col, TLSxy, flag = adjrowcol(df_pos, TLSxy)
        TLSedgexy = TLSxy.loc[TLSxy.TLS == sub_tls]
        
        bin_f = np.full((spot_col+2, spot_row+2), 0, dtype="uint8")
        bin_f[TLSedgexy['col'].values, TLSedgexy['row'].values] = 1
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))       
        erode = cv2.erode(bin_f, kernel)    
        edge = bin_f - erode               
        edgeXY = np.where(edge == 1)      

        Sizef = 5
        backcolor = "#eaeaea"       
        cm = plt.cm.Spectral_r

        fig, ax = plt.subplots(figsize=(Sizef * spot_row / spot_col, Sizef))
        ax.set_xlim((0, spot_row+3))
        ax.set_ylim((0, spot_col+3))
        r = 0.55
        r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
        marker_size = np.pi * r_**2
        for celli in range(len(Sel_CellTypeList)):
            sel_cell = Sel_CellTypeList[celli]
            subda = TLSedgexy.loc[TLSedgexy.CellSubType_new == sel_cell, ["row", "col"]]
            if len(subda) > 0:
                ax.scatter(x=subda['row']+1, y=subda['col']+1, c=mycolor_Plasma_only[celli], s=marker_size, edgecolors='black', linewidth=0, cmap=cm)

        plot_contour(TLSedgexy, spot_row, spot_col, 'slide', '#000000', 0.5)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        # plt.gca().add_patch(plt.Rectangle(xy=(0.5,0.5),
        #                     width=spot_row+2,
        #                     height=spot_col+2,
        #                     edgecolor='#bdbdbd',
        #                     fill=False, linewidth=1.25))
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0.1, wspace=0.1)
        ax.figure.savefig(f'{pathans}{sampid}/{sub_tls}_celltype_noFrame.pdf')
        print(sub_tls)