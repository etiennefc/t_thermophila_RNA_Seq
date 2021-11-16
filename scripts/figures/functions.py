#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.ticker
#import matplotlib.lines as lines
#import collections as coll
#from scipy import stats, misc

"""Functions to create various graphs"""

def heatmap_fixed(df, cmap, cbar_label, xtick_labels, path, **kwargs):
    """
    Creates a clustered heatmap horizontally but not vertically.
    """

    plt.rcParams['svg.fonttype'] = 'none'

    # Create the heatmap (clustered by rows and columns)
    graph = sns.clustermap(df, cmap=cmap, xticklabels=True, yticklabels=True, **kwargs)
    plt.xlabel(xlabel=cbar_label, fontsize=17)
    graph.ax_heatmap.xaxis.set_ticks_position('top')
    graph.ax_heatmap.set_xticklabels(xtick_labels)
    plt.setp(graph.ax_heatmap.xaxis.get_ticklabels(), fontsize=17)
    plt.setp(graph.ax_heatmap.yaxis.get_ticklabels(), fontsize=9)
    plt.savefig(path, bbox_inches='tight', dpi=900)


def pie_multiple(count_list, labels, colors, ax_title, title, path, **kwargs):
    """
    Create a pie chart for each sample type.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, axes = plt.subplots(1, 3, figsize=(25, 14))
    plt.subplots_adjust(hspace=0.5)
    ax = axes.flatten()
    for i, tissu in enumerate(count_list):
        count_per_tissue = count_list[i][:]
        ax[i].set_title(ax_title[i], fontdict={'fontsize': 25}, x=0.5, y=0.75)
        ax[i].pie(count_per_tissue, colors=colors, textprops={'fontsize': 25},
                    autopct=autopct_generator(0.5), **kwargs)
        ax[i].axis('equal')

    fig.suptitle(title, x=0.5, y=0.8, fontsize=30)
    fig.legend(labels=labels, loc='upper right', bbox_to_anchor=(0.6, 0.3), prop={'size': 30})
    plt.savefig(path, dpi=600)


def volcano(df, x_col, y_col, hue_col, xlabel, ylabel, title, color_dict, path,
            **kwargs):
    """
    Creates a violin plot (using a x, y and hue column).
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1,1, figsize=(10,10))

    #ax.set_xscale('symlog')
    #ax.set_yscale('symlog')
    #ax.set_xlim(-4, 8)
    #ax.set_ylim(-0.1, 250)

    #Set minor y tickmarks (linear between 0 and 1; log  with values > 1)
    #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * .1,
    #                                      numticks=100)
    #ax.yaxis.set_minor_locator(locmin)
    #ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    ax.set_title(title, fontsize=25)

    sns.scatterplot(data=df, x=x_col, y=y_col, hue=hue_col,
                    palette=color_dict, edgecolor='face',
                    s=25, **kwargs)

    plt.legend(prop={'size': 15})
    ax.set_xlabel(xlabel, fontsize=25)
    ax.set_ylabel(ylabel, fontsize=25)
    plt.savefig(path, bbox_inches='tight', dpi=600)


def scatter(df, x_col, y_col, hue_col, xlabel, ylabel, title, color_dict, path,
            **kwargs):
    """
    Creates a scatter plot (using a x, y and hue column).
    """

    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1,1, figsize=(10,10))
    ax.set_xscale('symlog')
    ax.set_yscale('symlog')
    ax.set_xlim(-10000000, 10000000)
    ax.set_ylim(-10000000, 10000000)
    ax.set_title(title, fontsize=25)
    #ax.axhline(y=0, color='grey', linestyle='--')
    #ax.axvline(x=0, color='grey', linestyle='--')
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color='grey', linestyle='--')  # add x=y diagonal
    sns.scatterplot(data=df, x=x_col, y=y_col, hue=hue_col,
                    palette=color_dict, edgecolor='face',
                    s=25, **kwargs)

    #plt.legend(prop={'size': 15})
    ax.set_xlabel(xlabel, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    plt.savefig(path, bbox_inches='tight', dpi=600)


def density_x(df_list, xlabel, ylabel, title, colors, crit_list, path, **kwargs):
    """
    Creates a density plot with a number x of dfs to represent on the same ax.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    for i, df in enumerate(df_list):
        sns.kdeplot(df, shade=True, ax=ax, color=colors[i], **kwargs)
    ax.set_xlabel(xlabel, fontdict={'fontsize': 25})
    ax.set_ylabel(ylabel, fontdict={'fontsize': 25})
    #ax.set_xscale('symlog')
    ax.set_xlim(-200, 60000)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)

    legend_list = []
    for i, crit in enumerate(crit_list):
        legend_element = mpatches.Patch(color=colors[i], label=crit)
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(1,1),
                fontsize=20)

    fig.suptitle(title, fontsize=20)
    plt.savefig(path, bbox_inches='tight', dpi=500)


def donut_2(counts, labels, colors, title, legend_labels, legend_colors, path,
            **kwargs):
    """
    Creates a donut plot with two layers. Counts, labels and colors are nested
    lists for outer and inner layers attributes.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots()
    ax.axis('equal')
    outer_donut, _, _ = ax.pie(counts[0], radius=1.3, labels=labels[0],
                        colors=colors[0], autopct='%.1f%%', pctdistance=0.85, **kwargs)
    plt.setp(outer_donut, width=0.4, edgecolor='white')

    inner_donut, _, _ = ax.pie(counts[1], radius=1.3-0.4, labels=labels[1],
                    colors=colors[1], autopct='%.1f%%', pctdistance=0.8, **kwargs)
    plt.setp(inner_donut, width=0.4, edgecolor='white')

    legend_list = []
    for i, crit in enumerate(legend_labels):
        legend_element = mpatches.Patch(color=legend_colors[i], label=crit)
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper left', bbox_to_anchor=(-0.1, 1),
                fontsize=10)

    plt.savefig(path, dpi=600)


def grouped_stacked_bar(lists, x_tick_labels, labels, title, xlabel,
                        ylabel, colors, legend_title, path, **kwargs):
    """
    Create a grouped stacked bar chart. Lists is a list of two lists of lists;
    labels is the labels of the stacked variables.
    """
    plt.rcParams['svg.fonttype'] = 'none'

    # Create dfs from the lists in 'lists'
    df1 = pd.DataFrame(lists[0], index=x_tick_labels, columns=labels)
    df2 = pd.DataFrame(lists[1], index=x_tick_labels, columns=labels)

    # Create the bar plot
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    df1.plot.bar(ax=ax, position=1, width=.3, color=colors, stacked=True,
                edgecolor='black', **kwargs)
    df2.plot.bar(ax=ax, position=0, width=.3, color=colors, stacked=True,
                edgecolor='black', **kwargs)
    plt.autoscale()
    plt.title(title)
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)


    # Add legend
    legend_list = []
    for i, crit in enumerate(labels):
        legend_element = mpatches.Patch(color=colors[crit], label=crit)
        legend_list.append(legend_element)
    legend = ax.legend(handles=legend_list, bbox_to_anchor=(1.1,1.1), fontsize=15)
    legend.set_title(legend_title,prop={'size':18})
    ax.add_artist(legend)

    plt.savefig(path, bbox_inches='tight', dpi=600)


def stacked_bar(lists, x_tick_labels, labels, title, xlabel, ylabel, colors, path, **kwargs):
    """
    Create a stacked bar chart from a list of lists ('lists').
    """
    rc = {'ytick.labelsize': 30, 'xtick.labelsize': 35}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'

    df = pd.DataFrame(lists, index=x_tick_labels, columns=labels)
    print(df)

    ax = df.plot.bar(stacked=True, figsize=(12,8), color=colors, **kwargs)
    ax.set_xticklabels(x_tick_labels, rotation=0)
    plt.legend(fontsize=35, loc=5, bbox_to_anchor=(0.5, 1.25))
    plt.title(title)
    plt.xlabel(xlabel, fontsize=40)
    plt.ylabel(ylabel, fontsize=40)
    plt.autoscale()
    plt.savefig(path, bbox_inches='tight', dpi=600)

def tpm_list(df, global_list, partial_col_title, crit_list, column_crit):
    """
    Creates a list of lists of abundance values (in TPM).
    """

    tpm_list = []
    for i, avg_sample in enumerate(global_list):
        temp_df_list = []
        for j, crit in enumerate(crit_list):
            temp_df = df[df[column_crit] == crit]
            total_tpm = temp_df[avg_sample+partial_col_title].sum()
            #avg_total_tpm = temp_df[tissu+titre_partiel_colonne].mean()
            temp_df_list.append(total_tpm)
            #temp_df_list.append(avg_total_tpm)
        tpm_list.append(temp_df_list)
    '''

    tpm_list = []
    for i, crit in enumerate(critère_list):
        temp_df_list = []
        for j, tissue in enumerate(tissu_liste):
            temp_df = df[df[colonne_crit] == crit]
            total_tpm = temp_df[tissue+titre_partiel_colonne].sum()
            temp_df_list.append(total_tpm)
        tpm_list.append(temp_df_list)
        '''
    return tpm_list



def count_list_x(initial_df, global_col, criteria, specific_col):
    """
    Create a list of lists using initial_col to split the global list and
    specific_col to create the nested lists.
    """
    df_list = []

    #Sort in acending order the unique values in global_col and create a list of
    # df based on these values
    print(sorted(list(initial_df[global_col].unique())))
    for val in sorted(list(initial_df[global_col].unique())):
        temp_val = initial_df[initial_df[global_col] == val]
        df_list.append(temp_val)


    l = []

    for i, df in enumerate(df_list):
        temp = []
        for j, temp1 in enumerate(criteria):
            crit = df[df[specific_col] == temp1]
            crit = len(crit)
            temp.append(crit)
        l.append(temp)

    return l


def percent_count(count_list):
    """
    Create from a list of lists a percentage list of lists.
    """
    percent_list = []

    for i, cat in enumerate(count_list):
        temp = []
        for j, crit in enumerate(cat):
            total = sum(cat)
            if total != 0:
                percent = crit/total * 100
                percent = round(percent, 2)
            else:
                percent = 0
            temp.append(percent)
        percent_list.append(temp)

    return percent_list

def autopct_generator(limit):
    """
    Show percentages (above 'limit') of all categories on a pie chart.
    """
    def inner_autopct(pct):
        return ('%1.0f%%' % pct) if pct > limit else ''
    return inner_autopct
