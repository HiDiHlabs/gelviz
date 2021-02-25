import matplotlib.pyplot as plt
import pybedtools
import pandas as pnd
import numpy as np
import tabix
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle
from matplotlib.patches import Arrow
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.cm as cm
import matplotlib
import tabix
import math

def plotGenes(genes_bed,
              exons_bed,
              introns_bed,
              region_bed,
              blacklist=None,
              gene_map=None,
              plot_gene_ids=True,
              y_max=None,
              distance_ratio=0.1,
              ax=None,
              plot_legend=False,
              legend_loc="lower right",
              color_plus="#80b1d3",
              color_minus="#fb8072"):
    """Function for plotting gene structures, i.e. introns exons of genes.

    :param genes_bed: :class:`pybedtools.BedTool` object containing TX start,
        and TX end of genes.
    :type genes_bed: :class:`pybedtools.BedTool`
    :param exons_bed: :class:`pybedtools.BedTool` object containing exons of
        genes.
    :type exons_bed: :class:`pybedtools.BedTool`
    :param introns_bed: :class:`pybedtools.BedTool` object containing introns
    :type introns_bed: :class:`pybedtools.BedTool`
    :param region_bed: :class:`pybedtools.BedTool` object containing the one
        region, for which the gene plot is created.
    :type region_bed: :class:`pybedtools.BedTool`
    :param blacklist: List of gene names, for genes that should not be shown on
        the plot, default is None
    :type blacklist: list, optional
    :param plot_gene_ids: If True, all gene ids will be included in the plot,
        False otherwise, default is True
    :type plot_gene_ids: bool, optional
    :param y_max: Max y value in the gene plot. If not set, then y_max is the
        max number of stacked genes, default is None.
    :type y_max: bool, optional
    :param distance_ratio: Minimal distance between two genes, as ratio of ax
        width, such that two genes are plotted side by side. If this ratio is
        underwent, the genes will be stacked, default is 0.1.
    :type distance_ratio: float, optional
    :param ax: Axes instance on which the genes are plotted, default is None.
    :type ax: :class:`matplotlib.axes._subplots.AxesSubplot`, optional
    :param plot_legend: If True, a legend describing plus or minus stranded
        genes is plotted, False otherwise. Default is False.
    :type plot_legend: bool, optional
    :param legend_loc: Location of the legend. Either of "lower left",
        "lower right", "upper left", "upper right", default is "lower right".
    :type legend_loc: str, optional
    :param color_plus: Color code for plus stranded genes, default is "#80b1d3".
    :type color_plus: str, optional.
    :param color_minus: Color code for minus stranded genes, default is
        "#fb8072".
    :type color_minus: str, optional.

    :return: Tuple of max_y_pos+1.5, patch_list, patch_description_list, where

        1. max_y_pos+1.5 is the max_y_position + 1.5. max_y_pos defines the \
number of stacked genes.
        2. patch_list is the list of patches drawn on the ax.
        3. patch_description_list is the list of descriptions for the patches \
drawn on the ax.

    :rtype: list
    """


    ax = ax if ax is not None else plt.gca()

    genes_in_region = genes_bed
    exons_in_region = exons_bed
    introns_in_region = introns_bed

    region_border_up = int(region_bed[0][1])
    region_border_down = int(region_bed[0][2])
    region_size = region_border_down-region_border_up

    color_forward = color_plus
    color_reverse = color_minus

    max_y_pos = None

    if(not len(genes_in_region) == 0):
        # Determine y positions of genes for plotting
        max_y_pos, y_pos_dict = determineYPosGene(genes_in_region,
                                                  (region_border_down-
                                                    region_border_up),
                                                  distance_ratio)

        if(not y_max is None):
            max_y_pos = y_max

        # Plot Exons
        for i in exons_in_region:
            start = int(i[1])
            end = int(i[2])
            gene_name = str(i[3])
            if(not blacklist is None and gene_map[gene_name] in blacklist):
                continue

            # Define color for gene plotting
            strand = str(i[5])
            color = color_forward
            if(strand == "-"):
                color = color_reverse
            y = max_y_pos-y_pos_dict[gene_name]+0.5

            rect = Rectangle((start, y-.2),
                             end-start,
                             .4,
                             color=color,
                             capstyle='butt',
                             linewidth=0)

            ax.add_patch(rect)

        patch_list = []
        patch_description_list = []
        met_forward = False
        met_reverse = False
        # Plot Introns
        for i in introns_in_region:
            start = int(i[1])
            end = int(i[2])
            gene_name = str(i[3])

            if(not blacklist is None and gene_map[gene_name] in blacklist):
                continue
            # Define color for gene plotting
            strand = str(i[5])
            color = color_forward
            if(strand == "-"):
                color = color_reverse
            y = max_y_pos-y_pos_dict[gene_name]+0.5
            patch = Rectangle((start, y-.03),
                              end-start,
                              .06,
                              color=color,
                              capstyle='butt',
                              linewidth=0)
            ax.add_patch(patch)

            if(strand == "+" and not(met_forward)):
                patch_list += [patch]
                patch_description_list += ["forward strand"]
                met_forward = True
            elif(strand == "-" and not(met_reverse)):
                patch_list += [patch]
                patch_description_list += ["reverse strand"]
                met_reverse = True

        # Plot Gene Names
        if(plot_gene_ids):
            for i in genes_in_region:
                start = int(i[1])
                gene_name = str(i[3])
                if(not blacklist is None and gene_map[gene_name] in blacklist):
                    continue

                # Define color for gene plotting
                strand = str(i[5])
                color = color_forward
                if(strand == "-"):
                    color = color_reverse

                border_distance_down = region_border_down-start
                if(start < region_border_up):
                    start = region_border_up
                    border_distance_down = region_border_down-start
                    if(not(float(border_distance_down)/float(region_size)
                           < distance_ratio)):
                        gene_name = str(i[3])
                        gene_name_label = gene_name
                        if(not gene_map is None):
                            gene_name_label = gene_map[gene_name]

                        y = max_y_pos-y_pos_dict[gene_name]+.8
                        plt.text(start,
                                 y,
                                 gene_name_label,
                                 size=5,
                                 color = color)
                gene_name = str(i[3])
                gene_name_label = gene_name
                if(not gene_map is None):
                    gene_name_label = gene_map[gene_name]

                y = max_y_pos-y_pos_dict[gene_name]+.8
                plt.text(start, y, gene_name_label, size=5, color = color)

        plt.xlim([region_border_up, region_border_down])
        plt.ylim([0, max_y_pos+1.5])
        plt.yticks([], [])

        if(plot_legend):
            plt.legend(patch_list,
                       patch_description_list,
                       loc=legend_loc,
                       fontsize=5)

        return max_y_pos+1.5, patch_list, patch_description_list

def determineYPosGene(genes_bed,
                      region_size,
                      distance_ratio):
    '''Function that determines the max y position for gene plotting via
    function plotGenes.

    :param genes_bed: :class:`pybedtools.BedTool` object containing genes to be
        plotted.
    :type genes_bed: :class:`pybedtools.BedTool`
    :param region_size: Size of region to be plotted in base pairs.
    :type region_size: int
    :param distance_ratio: Minimal distance between two genes, as ratio of ax
        width, such that two genes are plotted side by side. If this ratio is
        underwent, the genes will be stacked.
    :type distance_ratio: float

    :return: Tuple of

        1. max_y_pos: Defines the number of stacked genes.
        2. y_pos_dict: Dictionary with keys = gene ids and values = y position \
of gene.

    :rtype: tuple

    '''
    sort_indices = [int(idx) for idx in np.argsort([i[1] for i in genes_bed])]
    genes_sorted_bed = [genes_bed[i] for i in sort_indices]

    y_pos_dict = {}

    y_level_dict = {}

    max_y_pos = 0
    for interval in genes_sorted_bed:
        gene_name = interval[3]
        gene_start = int(interval[1])
        gene_end = int(interval[2])

        for i in range(max_y_pos+1):
            if(i == 0 and not max_y_pos in y_level_dict):
                y_pos_dict[gene_name] = i
                y_level_dict[i] = [[gene_start, gene_end]]
                break
            elif(gene_start > y_level_dict[i][-1][1] and
                 float(gene_start-y_level_dict[i][-1][0])/float(region_size) >
                 distance_ratio):
                y_pos_dict[gene_name] = i
                y_level_dict[i] += [[gene_start, gene_end]]
                break
            elif(i == max_y_pos):
                max_y_pos += 1
                y_pos_dict[gene_name] = max_y_pos
                y_level_dict[max_y_pos] = [[gene_start, gene_end]]
                break
            else:
                continue

    return max_y_pos, y_pos_dict

def createGeneNameMap(gene_name_mapping_filename):
    '''Function that creates a mapping between gene ids

    :param gene_name_mapping_file: Path to a tab separated file, for which the
        first column is a ensemble gene id, and the second column is the HUGO
        gene name
    :type gene_name_mapping_file: str

    :return: Dictionary containing the gene id mapping.
    :rtype: dictionary
    '''
    gene_name_mapping_file = open(gene_name_mapping_filename, "r")

    gene_map = {}

    for line in gene_name_mapping_file:
        split_line = line.rstrip().split("\t")

        ensembl_gene_id = split_line[0].split(".")[0]

        hugo_gene_symbol = split_line[1].split(".")[0]

        gene_map[ensembl_gene_id] = hugo_gene_symbol

    gene_name_mapping_file.close()

    return gene_map

def plotGeneExpression(genes_bed,
                       region_bed,
                       expression_df_g1,
                       expression_df_g2,
                       gene_names_map,
                       blacklist=None,
                       ax=None,
                       plot_legend=False,
                       color_g1="#fb8072",
                       color_g2="#80b1d3",
                       g1_id="tumor",
                       g2_id="normal",
                       plot_gene_names=True):
    '''Function for plotting paired gene expression (e.g. tumor and normal) on a
    gene region scale retaining the position of genes.

    :param genes_bed: :class:`pybedtools.BedTool` object containing TXstart,
        and TXend of genes.
    :type genes_bed: :class:`pybedtools.BedTool`
    :param region_bed: :class:`pybedtools.BedTool` object containing the region
        to be plotted
    :type region_bed: :class:`pybedtools.BedTool`
    :param expression_df_g1: :class:`pandas.Dataframe` containing the expression
        values of g1 samples (columns: sample ids; index: gene ids)
    :type expression_df_g1: :class:`pandas.DataFrame`
    :param expression_df_g2: :class:`pandas.Dataframe` containing the expression
        values of g2 samples (columns: sample ids; index: gene ids)
    :type expression_df_g2: :class:`pandas.DataFrame`
    :param gene_names_map: Dictionary with keys: ENSEMBL GENE IDs, and values:
        HUGO GENE SYMBOLs.
    :type gene_names_map: dict.
    :param blacklist: Set containing gene ids not to be plotted, default to
        None.
    :type blacklist: set, optional
    :param ax: Axis used for plotting, defaults to None.
    :type ax: :class:`matplotlib.axes._subplots.AxesSubplot`, optional
    :param plot_legend: If True legend is plotted, False otherwise, defaults to
        False.
    :type plot_legend: bool
    :param color_g1: Color used for plotting g1 samples expression, defaults to
        "#fb8072".
    :type color_g1: str, optional
    :param color_g2: Color used for plotting g2 samples expression, defaults to
        "#80b1d3".
    :type color_g2: str, optional
    :param g1_id: ID of g1 used for legend plotting, defaults to "tumor".
    :type g1_id: str, optional
    :param g2_id: ID of g2 used for legend plotting, defaults to "normal".
    :type g2_id: str, optional
    :param plot_gene_names: If True, the HUGO GENE SYMBOLs will be shown, else
        the GENE SYMBOLs are hidden.
    :type plot_gene_names: bool.

    :return: Axis on which plot was placed.
    :rtype: :class:`matplotlib.axes._subplots.AxesSubplot`
    '''
    ax = ax if ax is not None else plt.gca()

    # Get gene names and regions
    genes_in_region_bed = genes_bed.intersect(region_bed,
                                              wa=True,
                                              u=True).sort()

    gene_names = []
    gene_regions = []

    for e in genes_in_region_bed:
        gene_name_ens = str(e[3])
        gene_names += [gene_names_map[gene_name_ens]]
        gene_regions += [[int(e[1]), int(e[2])]]

    region_right_border = int(region_bed[0][2])
    region_left_border = int(region_bed[0][1])

    # Determine minimal extension of barplot
    extension=None
    for i in range(len(gene_regions)):
        if(not blacklist is None and gene_names[i] in blacklist):
            continue

        left_border = gene_regions[i][0]
        right_border = None
        if(i < len(gene_names)-1):
            right_border = gene_regions[i+1][0]
        else:
            right_border = region_right_border
        current_extension = right_border-left_border

        if(current_extension == 0.):
            continue

        if(extension is None):
            extension = float(current_extension)
        elif(current_extension < extension):
            extension = float(current_extension)

    boxprops = {"color": "k", "linewidth": .3}
    flierprops = {"color": "k"}
    medianprops = {"color": "k", "linewidth": .3}
    whiskerprops = {"color": "k", "linewidth": .3}
    capprops={"color": "k", "linewidth": .3}

    patch_list = None
    patch_description_list = None

    tick_positions = []
    gene_names_clean = []
    counter=0
    patch_saved = False
    for gene_name in gene_names:
        left_border = gene_regions[counter][0]
        right_border = region_right_border

        if(not blacklist is None and gene_name in blacklist):
            counter += 1
            continue

        if(counter < len(gene_names)-1):
            right_border = gene_regions[counter+1][0]

        bplot_g1_pos = left_border + extension/4.
        bplot_g2_pos = left_border + 3*(extension/4.)

        tick_positions += [left_border + extension/2.]
        gene_names_clean += [gene_name]
        exp_values_g1 = expression_df_g1.loc[gene_name, :]
        if(type(exp_values_g1).__name__ == "Series"):
            exp_values_g1 = list(exp_values_g1)
        else:
            exp_values_g1 = list(exp_values_g1.iloc[0, :])
        exp_values_g2 = expression_df_g2.loc[gene_name, :]
        if(type(exp_values_g2).__name__ == "Series"):
            exp_values_g2 = list(exp_values_g2)
        else:
            exp_values_g2 = list(exp_values_g2.iloc[0, :])

        bplot_g1 = ax.boxplot([np.log2([i if
                                        i >= 1. else
                                        1. for
                                        i in exp_values_g1])],
                              positions=[bplot_g1_pos],
                              widths=extension/2.,
                              patch_artist=True,
                              boxprops=boxprops,
                              flierprops=flierprops,
                              medianprops=medianprops,
                              whiskerprops=whiskerprops,
                              capprops=capprops,
                              showfliers=False)
        bplot_g2 = ax.boxplot([np.log2([i if
                                        i >= 1. else
                                        1. for
                                        i in exp_values_g2])],
                              positions=[bplot_g2_pos],
                              widths=extension/2.,
                              patch_artist = True,
                              boxprops=boxprops,
                              flierprops=flierprops,
                              medianprops=medianprops,
                              whiskerprops=whiskerprops,
                              capprops=capprops,
                              showfliers=False)

        bplot_g1["boxes"][0].set_facecolor(color_g1)
        bplot_g2["boxes"][0].set_facecolor(color_g2)

        if(not patch_saved):
            patch_saved=True
            patch_list = [bplot_g1["boxes"][0], bplot_g2["boxes"][0]]
            patch_description_list = [g1_id, g2_id]

        counter += 1

    ax.set_xlim(region_left_border, region_right_border)
    ax.xaxis.set_major_locator(ticker.FixedLocator((tick_positions)))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter((gene_names_clean)))
    if(not plot_gene_names):
        ax.xaxis.set_major_formatter(
            ticker.FixedFormatter(([ " " for i in
                                    gene_names_clean])))
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
        tick.set_size(6)
    for ytick in ax.get_yticklabels():
        ytick.set_size(6)

    if(plot_legend):
        ax.legend(patch_list,
                  patch_description_list,
                  fontsize=5,
                  loc='lower left')

    return ax

def plotGeneExpressionEqualDist(genes_bed,
                                gene_mid_points,
                                region,
                                expression_df,
                                groups,
                                gene_names_map=None,
                                blacklist=None,
                                ax=None,
                                plot_legend=False,
                                colors=None,
                                ids=None,
                                plot_gene_names=True,
                                position_gene_names="bottom",
                                log_transformed=True,
                                plot_points=False,
                                alpha=.5):
    '''Function for plotting grouped gene expression (e.g. tumor and normal) on
    a gene region scale equalizing the position of genes.

    :param genes_bed: :class:`pybedtools.BedTool` object containing gene
        regions.
    :type genes_bed: :class:`pybedtools.BedTool`
    :param gene_mid_points: list of integer values containing center positions
        of genes.
    :type gene_mid_points: list
    :param region: List containing the region to be plotted
        ([<chrom>, <start>, <end>]).
    :type region: list
    :param groups: List of lists containing the IDs of the different groups.
    :type groups: list
    :param gene_names_map: Dictionary with keys: ENSEMBL GENE IDs, and values:
        HUGO GENE SYMBOLs.
    :type gene_names_map: dict.
    :param expression_df: class:`pandas.DataFrame` object containing the
        expression values of all samples (columns: sample ids; index: gene ids).
    :type expression_df: class:`pandas.DataFrame`
    :param blacklist: Set containing gene ids not to be plotted, defaults to
        None,
    :type blacklist: set, optional
    :param ax: (default: None) Axis used for plotting, defaults to None.
    :type ax: :class:`matplotlib.axes._subplots.AxesSubplot`, optional
    :param plot_legend: If True plot legend, False otherwise, defaults to False.
    :type plot_legend: bool, optional
    :param colors: List of colors used for plotting samples expression. The
        number of colors must be the same as the number of groups, defaults to
        None.
    :type colors: str, optional
    :param ids: IDs used for legend plotting, defaults to None. Number of ids
        must be the same as the number of groups.
    :type ids: list, optional.
    :param plot_gene_names: True if gene names shall be plotted,
        False otherwise, defaults to True.
    :type plot_gene_names: bool, optional
    :param position_gene_names: Either of "top", or "bottom", defaults to
        "bottom".
    :type position_gene_names: str, optional
    :param log_transformed: If True use log transformed values for plotting,
        non-transformed values otherwise.
    :type log_transformed: bool, optional
    :param plot_points: If True, a point per expression value is plotted in
        addition to the boxplot, no points are plotted otherwise, defaults to
        False.
    :type plot_points: bool, optional
    :param alpha: Alpha value for the background color of the boxplots boxes,
        defaults to 0.5.
    :type alpha: float, optional

    :return: Plots axis.
    :rtype: :class:`matplotlib.axes._subplots.AxesSubplot`
    '''
    standard_colors = ["#66c2a5",
                       "#fc8d62",
                       "#8da0cb",
                       "#ec87c2",
                       "#a6d854",
                       "#ffd92f",
                       "#e5c494",
                       "#bbbbbb"]

    ax = ax if ax is not None else plt.gca()

    region_bed = pybedtools.BedTool("\t".join([str(i) for i in region]),
                                    from_string=True)

    # Get gene names and regions
    genes_in_region_bed = genes_bed.intersect(region_bed,
                                              wa=True,
                                              u=True).sort()

    gene_names = []
    gene_regions = []

    for e in genes_in_region_bed:
        gene_name_ens = str(e[3])
        if(not gene_names_map is None):
            gene_names += [gene_names_map[gene_name_ens]]
        else:
            gene_names += [gene_name_ens]
        gene_regions += [[int(e[1]), int(e[2])]]

    region_right_border = int(region_bed[0][2])
    region_left_border = int(region_bed[0][1])

    # Determine minimal extension of barplot
    extension=None
    if(len(gene_mid_points) <= 1):
        extension=region[2]-region[1]
    else:
        extension=gene_mid_points[1]-gene_mid_points[0]
    # Subtract a small percentage of region size from extension
    extension=extension-(region[2]-region[1])*.01

    boxprops = {"color": "k", "linewidth": .3, "alpha":alpha}
    flierprops = {"color": "k"}
    medianprops = {"color": "k", "linewidth": .3}
    whiskerprops = {"color": "k", "linewidth": .3}
    capprops={"color": "k", "linewidth": .3}

    patch_list = []
    patch_description_list = []

    tick_positions = []
    gene_names_clean = []
    counter=0
    for gene_name in gene_names:
        left_border = gene_mid_points[counter]-extension/2
        right_border = gene_mid_points[counter]+extension/2

        if(not blacklist is None and gene_name in blacklist):
            counter += 1
            continue

        n_groups = len(groups)
        for g in range(n_groups):
            bplot_pos = left_border + (2*g+1)*extension/float((n_groups*2.))

            tick_positions += [left_border + extension/2.]
            gene_names_clean += [gene_name]
            exp_values = expression_df.loc[gene_name, groups[g]]
            if(type(exp_values).__name__ == "Series"):
                exp_values = list(exp_values)
            else:
                exp_values = list(exp_values.iloc[0, :])

            expression_values = exp_values
            if(log_transformed):
                expression_values = np.log2([i
                                             if i >= 1.
                                             else 1.
                                             for i in exp_values])
            bplot = ax.boxplot(expression_values,
                               positions=[bplot_pos],
                               widths=extension/float(n_groups),
                               patch_artist=True,
                               boxprops=boxprops,
                               flierprops=flierprops,
                               medianprops=medianprops,
                               whiskerprops=whiskerprops,
                               capprops=capprops,
                               showfliers=False)

            color = None
            if(not colors is None):
                color = colors[g]
            else:
                color = standard_colors[g]

            bplot["boxes"][0].set_facecolor(color)

            if(plot_points):
                x_positions = [ (bplot_pos+
                                 (i-.5)*
                                 ((2*extension)/(float(n_groups)*3))) for i in
                               list(np.random.rand(len(expression_values))) ]
                plt.plot(x_positions, expression_values, "k.", markersize=3)

            g_id = None
            if(not ids is None):
                g_id = ids[g]
            else:
                g_id = "group "+str(g)

            if(not g_id in patch_description_list):
                patch_list += [bplot["boxes"][0]]
                patch_description_list += [g_id]

        counter += 1

    ax.set_xlim(region_left_border, region_right_border)
    if(position_gene_names == "top"):
        ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_major_locator(ticker.FixedLocator((tick_positions)))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter((gene_names_clean)))
    if(not plot_gene_names):
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(
            ([ " " for i in
              gene_names_clean])))
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
        tick.set_size(5)
    for ytick in ax.get_yticklabels():
        ytick.set_size(5)

    if(plot_legend):
        ax.legend(patch_list,
                  patch_description_list,
                  fontsize=5,
                  loc='lower left')

    return ax

def plotGenomicSegments(segments_list,
                        chrom,
                        start,
                        end,
                        ax = None):
    '''Function for plotting genomix segments in different colors

    :param segments_tabix_filename: Path to tabixed bed file containing
        (chrom, start, end, name, score, strand, start, end, color). The color
        field is used to determine the color for plotting (R,G,B).
    :type segments_Tabix_filename: str
    :param chrom: Chromosome of the region to be plotted.
    :type chrom: str
    :param start: Start position of the region to be plotted.
    :type start: str
    :param end: End position of the region to be plotted.
    :type end: str
    :param ax: Axis used for plotting, defaults to None.
    :type ax: :class:`matplotlib.axes._subplots.AxesSubplot`, optional

    :return: Dictionary with keys = names of segments, and values patch
    :rtype: dict
    '''

    ax = ax if ax is not None else plt.gca()

    patches_dict = {}
    for segment in segments_list:
        segment_start = int(segment[1])
        segment_end = int(segment[2])
        color = tuple([ float(i)/256. for i in
                       str(segment[-1]).split(",") ]+[1])
        segment_type = str(segment[3])
        if(segment_type == "R"):
            color = (1,1,1,1)

        rect = Rectangle((segment_start, 0),
                         segment_end-segment_start,
                         1,
                         color=color)
        ax.add_patch(rect)

        patches_dict[segment_type] = rect

    plt.xlim(int(start), int(end))
    plt.ylim(0, 1)
    plt.yticks([], [])

    return patches_dict

def plotCNVs(cnvs_bed,
             chromosome,
             start,
             end,
             ploidy=2,
             cnv_threshold=0.7,
             color_gain="g",
             color_loss="r",
             color_neutral="k",
             ax=None):
    '''Function for plotting CNV segments

    :param cnvs_bed: :class:`pybedtools.BedTool` object containing CNVs with
        following entries:

        1. Chromosome,
        2. Start Position,
        3. End Position,
        4. Deviation from ploidy,
        5. True Copy Number)

    :type cnvs_bed: :class:`pybedtools.BedTool`
    :param chromosome: Chromosome for which to plot CNVs.
    :type chromosome: str
    :param start: Start position on chromosome.
    :type start: int
    :param end: End position on chromosome.
    :type end: int
    :param ploidy: Assumed ploidy of tumor, defaults to 2.
    :type ploidy: int, optional
    :param cnv_threshold: Minimal deviation from ploidy to be considered as a
        CNV, defaults to 0.7.
    :type cnv_threshold: float, optional
    :param color_gain: Plot color of copy number gains, defaults to "g".
    :type color_gain: str, optional
    :param color_loss: Plot color of copy number losses, defaults to "r".
    :type color_loss: str, optional
    :param color_neutral: Plot color of copy number neutral regions, defaults to
        "k".
    :type color_neutral: str, optional
    :param ax: Axis used for plotting.
    :type ax: :class:`matplotlib.axes._subplots.AxesSubplot`, optional

    :return: Nothing to be returned
    :rtype: None
    '''
    # Use given axis for plotting
    ax = ax if ax is not None else plt.gca()

    for interval in cnvs_bed:
        current_start = int(interval[1])
        current_end = int(interval[2])
        ploidy_dev = float(interval[3])
        tcn = float(interval[4])

        # Smooth tcn, if ploidy_dev is smaller than cnv_threshold
        if(abs(ploidy_dev) < cnv_threshold):
            tcn = ploidy

        color = color_neutral
        if(ploidy_dev >= cnv_threshold):
            color=color_gain
        elif(ploidy_dev <= -1.*cnv_threshold):
            color = color_loss

        if(abs(ploidy_dev) > cnv_threshold):
            rect = Rectangle((current_start, tcn-.2),
                             current_end-current_start,
                             .4,
                             color=color,
                             edgecolor='none',
                             capstyle='butt',
                             linewidth=0)
            ax.add_patch(rect)
        else:
            rect = Rectangle((current_start, tcn-.1),
                             current_end-current_start,
                             .2,
                             color=color,
                             edgecolor='none',
                             capstyle='butt',
                             linewidth=0)
            ax.add_patch(rect)

    # Plot thresholds
    color_threshold=(189./255., 189./255., 189./255., 0.5)
    if(ploidy == 2):
        plt.plot([int(start), int(end)],
                 [1, 1],
                 color=color_threshold,
                 linestyle="--",
                 linewidth=.5)
        plt.plot([int(start), int(end)],
                 [2, 2],
                 color=color_threshold,
                 linestyle="--",
                 linewidth=.5)
        plt.plot([int(start), int(end)],
                 [3, 3],
                 color=color_threshold,
                 linestyle="--",
                 linewidth=.5)
        plt.plot([int(start), int(end)],
                 [4, 4],
                 color=color_threshold,
                 linestyle="--",
                 linewidth=.5)
    elif(ploidy == 4):
        plt.plot([int(start), int(end)],
                 [1, 1],
                 color=color_threshold,
                 linestyle="--",
                 linewidth=.5)
        plt.plot([int(start), int(end)],
                 [2, 2],
                 color=color_threshold,
                 linestyle="--",
                 linewidth=.5)
        plt.plot([int(start), int(end)],
                 [3, 3],
                 color=color_threshold,
                 linestyle="--",
                 linewidth=.5)
        plt.plot([int(start), int(end)],
                 [4, 4],
                 color=color_threshold,
                 linestyle="--",
                 linewidth=.5)
        plt.plot([int(start), int(end)],
                 [5, 5],
                 color=color_threshold,
                 linestyle="--",
                 linewidth=.5)
        plt.plot([int(start), int(end)],
                 [6, 6],
                 color=color_threshold,
                 linestyle="--",
                 linewidth=.5)

    plt.xlim([int(start), int(end)])
    if(ploidy == 2):
        plt.ylim([0, 4.5])
        plt.yticks([0, 1, 2, 3, 4], ["0", "1", "2", "3", "4"], size=6)
    elif(ploidy == 4):
        plt.ylim([0, 6.5])
        plt.yticks([0, 2, 4, 6], ["0", "2", "4", "6"], size=6)
    plt.xticks(rotation=45)

def plotCNVsHeat(cnvs_bed,
                 chromosome,
                 start,
                 end,
                 ploidy=2,
                 cnv_threshold=0.7,
                 cmap="bwr",
                 max_dev=None,
                 ax=None):
    '''Function for plotting CNV segments as heatmap

    :param cnvs_bed: :class:`pybedtools.BedTool` object containing CNVs with
        following entries:

        1. Chromosome,
        2. Start Position,
        3. End Position,
        4. Deviation from ploidy,
        5. True Copy Number)

    :type cnvs_bed: :class:`pybedtools.BedTool`
    :param chromosome: Chromosome for which to plot CNVs.
    :type chromosome: str
    :param start: Start position on chromosome.
    :type start: int
    :param end: End position on chromosome.
    :type end: int
    :param ploidy: Assumed ploidy of tumor, defaults to 2.
    :type ploidy: int, optional
    :param cnv_threshold: Minimal deviation from ploidy to be considered as a
        CNV, defaults to 0.7.
    :type cnv_threshold: float, optional
    :param cmap: Colormap used for plotting CNVs, defaults to "bwr".
    :type cmap: str, optional
    :param max_dev: Maximal deviation from ploidy to plot, defaults to None.
    :type max_dev: float, optional
    :param ax: Axis used for plotting, defaults to None.
    :type ax: :class:`matplotlib.axes._subplots.AxesSubplot`, optional

    :return: Nothing to be returned.
    :rtype: None
    '''
    # Use given axis for plotting
    ax = ax if ax is not None else plt.gca()

    colors = plt.cm.get_cmap(cmap)

    if(max_dev is None):
        max_dev = max([abs(float(i[3])) for i in cnvs_bed])

    for interval in cnvs_bed:
        current_start = int(interval[1])
        current_end = int(interval[2])
        ploidy_dev = float(interval[3])
        tcn = float(interval[4])

        if(tcn < -1.*max_dev):
            tcn = -1.*max_dev
        elif(tcn > max_dev):
            tcn = max_dev

        color = colors((ploidy_dev+max_dev)/(2*max_dev))

        if(abs(ploidy_dev) < cnv_threshold):
            color=colors(.5)

        rect = Rectangle((current_start, .5),
                         current_end-current_start,
                         1,
                         color=color,
                         edgecolor='none',
                         capstyle='butt',
                         linewidth=0)
        ax.add_patch(rect)

    plt.xlim([int(start), int(end)])
    plt.ylim([.5, 1.5])
    plt.xticks([], [])
    plt.yticks([], [])

def readACESeqAsBed(input_filename):
    '''Function that reads CNVs from ACESeq ("*most_important*") files and
    converts them to pybedtools.BedTool object

    :param input_filename: Full path to ACESeq "most_important" file
    :type input_filename: str

    :return: :class:`pybedtools.BedTool` object containing CNVs from ACESeq
    :rtype: :class:`pybedtools.BedTool`
    '''
    input_file = open(input_filename, "r")

    cnv_bed_list = []

    ploidy = None
    for line in input_file:
        if(line[:7] == "#ploidy"):
            ploidy = float(line.rstrip().split(":")[1])
            print(ploidy)

        if(line[0] == "#" or line[:5] == "chrom"):
            continue

        split_line = line.rstrip().split("\t")

        ploidy_dev = float(split_line[5])-ploidy

        chrom = split_line[0]
        if(chrom == "23"):
            chrom="X"
        elif(chrom == "24"):
            chrom = "Y"

        cnv_bed_list += [ [chrom,
                           split_line[1],
                           split_line[2],
                           str(ploidy_dev),
                           split_line[5],
                           "+"]
                        ]

    input_file.close()

    return pybedtools.BedTool("\n".join(["\t".join(e) for e in
                                         cnv_bed_list]),
                              from_string=True)

def plotChIPSignals(chip_signals,
                    r_chrom,
                    r_start,
                    r_end,
                    ax=None,
                    color="b",
                    offset=None,
                    merge=None):
    '''Function that plots bedGraph like iterators.

    :param chip_signals: Iterator for which each element is a list-ike
        object containing:

        1. Chromosome
        2. Start postion
        3. End position
        4. Value to be plotted as bar

    :type chip_signals: iterator
    :param r_chrom: Chromosome of region to be plotted.
    :type r_chrom: str
    :param r_start: Start position of region to be plotted.
    :type r_start: int
    :param r_end: End position of region to be plotted.
    :type r_end: int
    :param ax: Axis of plot
    :type ax: :class:`matplotlib.axes._subplots.AxesSubplot`, optional
    :param color: color of bars, defaults to "b".
    :type color: str, optional
    :param offset: Length of intervals, defaults to None.
    :type offset: int, optional
    :param merge: Number of elements to be merged. If this value is not equal to
        0, than merge elements will be averaged an plotted, defaults to 0.
    :type merge: int, optional

    :return: Nothing to be returned.
    :rtype: None

    '''

    ax = ax if ax is not None else plt.gca()

    max_signal = 0
    left = []
    height = []
    for signal in chip_signals:
        start = int(signal[1])
        end = int(signal[2])
        value = float(signal[3])
        if(value > max_signal):
            max_signal = value

        if(not offset is None):
            end = start + offset

        left += [start]
        height += [value]

    left_merged = []
    height_merged = []
    if(not merge is None):
        heights = []
        lefts = []
        for i in range(len(left)):
            if(i % merge == 0 and not (i == 0)):
                left_merged += [lefts[0]]
                lefts = []
                height_merged += [np.mean(heights)]
                heights = []
            heights += [height[i]]
            lefts += [left[i]]

        if(not i % merge == 0):
            left_merged += [lefts[0]]
            lefts = []
            height_merged += [np.mean(heights)]
            heights = []
        offset = merge*offset
        left = left_merged
        height = height_merged

    plt.bar(left, height, offset, color = color, edgecolor = color)

    plt.xlim(r_start, r_end)

def plotMethylationProfileHeat(methylation_bed,
                               chrom,
                               start,
                               end,
                               bin_size=1000,
                               ax = None):
    '''Function for plotting methylation values as heatmap

    :param methylation_bed: Methylation calls. Following fields must be
        included: Chrom, Start, End, Methylated Cs, Unmethylated Cs.
    :type methylation_bed: :class:`pybedtools.BedTool`
    :param chrom: Chromosome of region to be plotted.
    :type chrom: str
    :param start: Start position of region to be plotted.
    :type start: int
    :param end: End position of region to be plotted.
    :type end: int
    :param bin_size: size of bin to average methylation values, defaults to
        1000.
    :type bin_size: int, optional
    :param ax: Axis to be used for plotting, defaults to None.
    :type ax: :class:`matplotlib.axes._subplots.AxesSubplot`, optional

    :return: Nothing to be returned
    :rtype: None
    '''

    ax = ax if ax is not None else plt.gca()

    binned_meth_calls = [ [0, 0] for i in range(int(((end-start)/bin_size)+1)) ]

    counter = 0
    for element in methylation_bed:
        # Determine bin
        position = int(element[1])
        if(position < start or position > end):
            continue
        n_meth = int(element[3])
        n_unmeth = int(element[4])
        current_bin = int((position-start)/bin_size)

        counter += 1

        binned_meth_calls[current_bin][0] += n_meth
        binned_meth_calls[current_bin][1] += n_unmeth

    binned_average_meth = [ float(i[0])/(float(i[0])+float(i[1]))
                            if (float(i[0])+float(i[1])) > 0
                            else "NA"
                            for i in binned_meth_calls ]

    binned_average_meth_no_missing = []
    n = len(binned_average_meth)
    for i in range(n):
        if(not binned_average_meth[i] == "NA"):
            binned_average_meth_no_missing += [binned_average_meth[i]]
        else:
            meth_before = (binned_average_meth[i-1]
                           if not i == 0
                           else "NA")
            meth_after = (binned_average_meth[i+1]
                          if not i == len(binned_average_meth)-1
                          else "NA")
            average_list = [ j
                             for j
                             in [meth_before, meth_after]
                             if not j == "NA" ]
            binned_average_meth_no_missing += [ (float(sum(average_list))/
                                                float(len(average_list)))
                                                if len(average_list) > 0
                                                else 0. ]

    binned_average_meth = binned_average_meth_no_missing

    # Plot average methylation values per bin
    # Define Colormap
    cmap = cm.bwr
    norm = matplotlib.colors.Normalize(vmin=0., vmax=1.)
    m = matplotlib.cm.ScalarMappable(norm = norm, cmap = cmap)
    for cbin in range(len(binned_average_meth)):
        rect = Rectangle((start+cbin*bin_size, 0),
                         bin_size,
                         1,
                         color=m.to_rgba(binned_average_meth[cbin]))
        ax.add_patch(rect)

    plt.xlim([start, end])
    plt.ylim([0, 1])
    plt.xticks([], [])
    plt.yticks([], [])

def plotMethylationProfile(meth_calls,
                           chrom,
                           start,
                           end,
                           color="k",
                           ax=None):
    '''Function that plots methylation values as dot plots.

    :param meth_calls: Iterator containing list-like elements with the following
        entries:

        1. Chromsome
        2. Start position
        3. end position
        4. Number methylated cytosines
        5. Number unmethylated cytosines

    :type meth_calles: iterator
    :param chrom: Chromosome of region to be plotted.
    :type chrom: str
    :param start: Start position of region to be plotted.
    :type start: int
    :param end: End position of region to be plotted.
    :type end: int
    :param color: Color of points representing methylation values, defaults to
        "k".
    :type color: str, optional
    :param ax: Axis of plot, defaults to None.
    :type ax: :class:`matplotlib.axes._subplots.AxesSubplot`, optional

    :return: Nothing to be returned
    :rtype: None

    '''
    ax = ax if ax is not None else plt.gca()

    plt.plot([ (float(m[1])+float(m[2]))/2. for m in meth_calls ],
             [ float(m[3])/(float(m[3])+float(m[4]))
               if not(float(m[3])+float(m[4]) == 0.)
               else 0. for m in meth_calls],
             color=color,
             marker=".",
             linestyle='None',
             markersize=1,
             alpha=.5)
    plt.ylim([0, 1])
    plt.xticks([], [])
    plt.xlim([start, end])

def plotTX(chrom_r,
           start_r,
           end_r,
           TX_pos,
           direction="right",
           color="k",
           ax=None):
    '''Function that plots a translocation event as a bar, showing the part
    of the genome that is translocated.

    :param chrom_r: Chromosome of the region to be plotted.
    :type chrom_r: str
    :param start_r: Start position of the region to be plotted.
    :type start_r: int
    :param end_r: End position of the region to be plotted.
    :type end_r: int
    :param TX_pos: Position of the translocation.
    :type TX_pos: int
    :param direction: Direction of the genomic part that is translocated. Either
        of "left" (upstream), or "right" (downstream), defaults to "left".
    :type direction: str, optional
    :param color: Color of the bar representing the translocation, defaults to
        "k".
    :type color: str, optional
    :param ax: Axis of plot, defaults to None.
    :type ax: :class:`matplotlib.axes._subplots.AxesSubplot`, optional

    :return: Nothing to be returned.
    :rtype: None
    '''
    ax = ax if ax is not None else plt.gca()

    TX_start = TX_pos
    TX_end = end_r
    if(direction == "left"):
        TX_start = start_r
        TX_end = TX_pos

    rect = Rectangle((TX_start, .4),
                     TX_end-TX_start,
                     .2,
                     color=color,
                     capstyle='butt',
                     linewidth=0)
    ax.add_patch(rect)

    plt.xlim([start_r, end_r])
    plt.ylim([0.3, 0.7])

def plotRegions(regions,
                start,
                end,
                color="#cbebc4",
                edgecolor=False,
                alpha=1,
                ax = None):
    '''Functions that plots genomic regions as simple rectangles.

    :param regions: Iterator containig list-like elements with the following
        entries:

        1. Chromosome
        2. Start position
        3. End position

    :type regions: iterator
    :param start: Start position of the region to be plotted.
    :type start: int
    :param end: End position of the region to be plotted.
    :type end: int
    :param color: Color of the rectangles representing the regions to be
        plotted, defaults to "#cbebc4".
    :type color: str, optional
    :param edge_color: Color of region edge. If False, no edge is plotted,
        defaults to False.
    :type edge_color: str, optional
    :param alpha: Alpha value of the rectangle, representing the region to be
        plotted, defaults to 1.
    :type alpha: float, optional.
    :param ax: Axis of plot, defaults to None.
    :type ax: :class:`matplotlib.axes._subplots.AxesSubplot`, optional

    :return: Nothing to be returned
    :rtype: None
    '''
    ax = ax if ax is not None else plt.gca()

    c = 0
    for region in regions:
        if(not edgecolor):
            current_color = color
            rect = Rectangle([int(region[1]), -.75],
                             int(region[2])-int(region[1]),
                             1.5,
                             facecolor=current_color,
                             edgecolor='none',
                             alpha=alpha)
            c += 1
        else:
            current_color = color
            rect = Rectangle([int(region[1]), -.75],
                             int(region[2])-int(region[1]),
                             1.5,
                             facecolor=current_color,
                             edgecolor=edgecolor,
                             alpha=alpha)
            c += 1

        ax.add_patch(rect)

    plt.xticks([], [])
    plt.yticks([], [])
    plt.xlim([start, end])
    plt.ylim([-1, 1])

def plotMotifDirections(motifs_bed,
                        start,
                        end,
                        head_width=0.2,
                        head_length=1000,
                        overhang=0,
                        color_plus="#80b1d3",
                        color_minus="#fb8072",
                        ax=None):
    '''Function that plots TF motifs as arrows, indicating their directionality.

    :param motifs_bed: :class:`pybedtools.BedTool` object containing regions
        of the TF sited to be plotted.
    :type motifs_bed: :class:`pybedtools.BedTool`
    :param start: Start position of the region to be plotted.
    :type start: int
    :param end: End position of the region to be plotted.
    :type end: int
    :param head_width: Width of the arrow head as proportion of the arrow,
        defaults to 0.2
    :type head_width: float, optional
    :param head_length: Length of the arrow in bp (depends on the region that
        is plotted), defaults to 1000.
    :type head_length: int, optional
    :param overhang: Fraction that the arrow is swept back (0 overhang means
        triangular shape). Can be negative or greater than one. Defaults to 0.
    :type overhang: float, optional
    :param color_plus: Color of plus stranded TF regions, defaults to "#80b1d3".
    :type color_plus: str, optional
    :param color_minus: Color of plus stranded TF regions, defaults to
        "#fb8072".
    :type color_minus: str, optional
    :param ax: Axis on which to plot contact map, defaults to None.
    :type ax: :class:`matplotlib.axes._subplots.AxesSubplot`, optional

    :return: Nothing to be returned.
    :rtype: None
    '''
    ax = ax if ax is not None else plt.gca()

    for motif in motifs_bed:
        motif_start = int(motif[1])
        motif_end = int(motif[2])
        strand = str(motif[3])

        arrow_start = motif_start
        arrow_end = motif_end
        color=color_plus
        dx = head_length
        if(strand == "-"):
            arrow_start = motif_end
            arrow_end = motif_start
            color = color_minus
            dx = -1.*head_length

        plt.arrow(arrow_start,
                  .5,
                  dx,
                  0,
                  head_width=head_width,
                  head_length=head_length,
                  overhang=overhang,
                  head_starts_at_zero=False,
                  edgecolor="none",
                  facecolor=color,
                  length_includes_head=True)

    plt.xlim([start, end])
    plt.ylim([0.4, 0.6])

def plotHiCContactMap(contact_map,
                      start,
                      end,
                      segment_size,
                      cmap="Greys",
                      vmin=None,
                      vmax=None,
                      location="top",
                      ax=None):
    '''Function that plots HiC contact maps as pyramid plots

    :param contact_map: Matrix that contains the intensity values of HiC
        contacts.
    :type contact_map: :class:`pandas.DataFrame`
    :param start: Chromosomal start position of region to be plotted.
    :type start: int
    :param end: Chromosomal end position of region to be plotted.
    :type end: int
    :param segment_size: Size of the segments for which contacts were called.
    :type segment_size: int
    :param cmap: Name of the colormap to be used for plotting HiC intensities,
        defaults to "Greys".
    :type cmap: str, optional
    :param vmin: Minimal value of intensity range to be plotted, defaults to
        None
    :type vmin: float, optional
    :param vmax: Maximal value of intensity range to be plotted, defaults to
        None.
    :type vmax: float, optional
    :param location: Either of "top" | "bottom". If location == "top", the
        pyramid points upwards, else if location == "bottom" the pyramid points
        downwards, defaults to top,
    :type location: str, optional
    :param ax: Axis on which to plot contact map, defaults to None.
    :type ax: :class:`matplotlib.axes._subplots.AxesSubplot`, optional

    :return: Nothing to be returned.
    :rtype: None
    '''

    ax = ax if ax is not None else plt.gca()

    contact_map_index1 = (start)/segment_size
    contact_map_index2 = ((end)/segment_size)+1

    sliced_contact_map = contact_map.iloc[contact_map_index1:contact_map_index2,
                                          contact_map_index1:contact_map_index2]

    if(vmin is None):
        vmin = 0
    if(vmax is None):
        vmax = np.percentile(contact_map, 99.9)

    colormap = plt.get_cmap(cmap)

    for i in range(contact_map_index1, contact_map_index2):
        y_range = (range(contact_map_index1+(i-contact_map_index1),
                        contact_map_index2)
                   if location == "top"
                   else range(contact_map_index1,
                              contact_map_index2-(contact_map_index2-i)))
        for j in y_range:
            # Define midpoint of rectangle
            midpoint = (i*segment_size+(j*segment_size-i*segment_size)/2.,
                        (j*segment_size-i*segment_size)/2.)
            vertices = [(midpoint[0]-segment_size/2., midpoint[1]),
                        (midpoint[0], midpoint[1]-segment_size/2.),
                        (midpoint[0]+segment_size/2., midpoint[1]),
                        (midpoint[0], midpoint[1]+segment_size/2.),
                        (midpoint[0]-segment_size/2., midpoint[1])
                        ]

            codes = [Path.MOVETO,
                    Path.LINETO,
                    Path.LINETO,
                    Path.LINETO,
                    Path.CLOSEPOLY,
                    ]

            path = Path(vertices, codes)
            intensity_value = contact_map.iloc[i, j]
            intensity_value = (intensity_value/vmax
                               if intensity_value <= vmax
                               else 1.)
            facecolor = colormap(intensity_value)
            patch = matplotlib.patches.PathPatch(path,
                                                 facecolor=facecolor,
                                                 edgecolor='none')
            ax.add_patch(patch)

    ax.set_xlim(start, end)
    if(location == "top"):
        ax.set_ylim(0, (end-start)/2.)
    else:
        ax.set_ylim(-1.*(end-start)/2., 0)

def distanceEqualizer(genomic_segments,
                      start,
                      end,
                      direction="top_down",
                      color="k",
                      ax = None):
    '''Function that plots arcs from unequal distances of genomic segments to
    equal distances.

    :param genomic_segments: List of segments for which distances shall be
        equalized (each segment is of the form [<chrom>, <start>, <end>])
    :type genomic_segments: list
    :param start: Start position of the genomic region.
    :type start: int
    :param end: End position of the genomic region.
    :type end: int
    :param color: Color of lines equalizing distances, defaults to "k".
    :type color: str, optional
    :param direction: Direction of distance equalization (top_down | bottom_up),
        defaults to "top_down".
    :type direction: str, optional.
    :param ax: Axis on which to plot, defaults to None.
    :type ax: :class:`matplotlib.axes._subplots.AxesSubplot`, optional

    :return: List of equalized region midpoints.
    :rtype: list
    '''

    ax = ax if ax is not None else plt.gca()

    # Calculate midpoints of original and distance equalized segments
    n_segments = len(genomic_segments)
    equalized_region_size = (end-start)
    if(n_segments > 0):
        equalized_region_size=(end-start)/n_segments

    equalized_region_mid_points = []
    for i in range(1, n_segments+1):
        equalized_region_mid_points += [((start+
                                         i*equalized_region_size)-
                                        equalized_region_size/2)]

    region_mid_points = []
    for e in genomic_segments:
        if(int(e[1]) < start):
            region_mid_points += [start+(int(e[2])-start)/2]
        elif(int(e[2]) > end):
            region_mid_points += [int(e[1])+(end-int(e[1]))/2]
        else:
            region_mid_points += [int(e[1])+(int(e[2])-int(e[1]))/2]

    for i in range(len(region_mid_points)):
        region_mid_point = region_mid_points[i]
        equalized_region_mid_point = equalized_region_mid_points[i]

        codes = []
        vertices = []

        if(direction == "top_down"):
            codes = [Path.MOVETO,
                    Path.LINETO,
                    Path.LINETO,
                    Path.LINETO]

            vertices = [(region_mid_point, 1),
                        (region_mid_point, .8),
                        (equalized_region_mid_point, .2),
                        (equalized_region_mid_point, 0)]
        else:
            codes = [Path.MOVETO,
                    Path.LINETO,
                    Path.LINETO,
                    Path.LINETO]

            vertices = [(region_mid_point, 0),
                        (region_mid_point, .2),
                        (equalized_region_mid_point, .8),
                        (equalized_region_mid_point, 1)]

        path = Path(vertices, codes)

        path_patch = PathPatch(path,
                               facecolor="none",
                               edgecolor=color,
                               linewidth=.5)

        ax.add_patch(path_patch)

    ax.axis("off")
    plt.xlim([start, end])
    plt.ylim([0, 1])

    return equalized_region_mid_points

def plotCoordinates(chrom,
                    start,
                    end,
                    color="k",
                    ax = None,
                    upper=True,
                    loc_coordinates="up",
                    revert_coordinates=False,
                    rotation=0):
    '''Function that plots genomic coordinates in a linea fashion.

    :param chrom: Chromosome of the region to be plotted.
    :type chrom: str
    :param start: Start position of the region to be plotted.
    :type start: int
    :param end: End position of the region to be plotted.
    :type end: int
    :param color: Color of the genomic scales elements, defaults to "k".
    :type color: str, optional
    :param ax: Axis of plot, defaults to None.
    :type ax: :class:`matplotlib.axes._subplots.AxesSubplot`, optional
    :param upper: If True, make less ticks, else if False make more ticks.
    :type upper: bool, optional
    :param loc_coordinates: Either of "up" | "down". If "up", plot ticks to
        upper direction, else if "down", plot ticks to lower direction, defaults
        to "up".
    :type loc_coordinates: str, optional
    :param revert_coordinates: If True, coordinates are reverted to decreasing
        order. Else, coordinates stay in increasing order, defaults to False.
    :type revert_coordinates: bool, optional
    :param rotation: Rotational angle of coordinate strings, defaults to 0.
    :type rotation: int, optional

    :return: Nothing to be returned.
    :rtype: None
    '''
    ax = ax if ax is not None else plt.gca()

    tick_size = 10**math.ceil((np.log10((end-start)/10)))
    if(not upper):
        tick_size = 10**int((np.log10((end-start)/10)))

    # Determine first tick position
    first_tick = start+(tick_size-start%tick_size)

    ticks = []

    current_tick = first_tick
    while(current_tick <= end):
       ticks += [current_tick]

       current_tick = current_tick + tick_size

    scale = None
    if(first_tick > 1000000):
        scale = "Mb"
    else:
        scale="Kb"

    digits_to_round = None
    divisor = None
    if(scale == "Mb"):
        digits_to_round = int(6-np.log10(tick_size))
        divisor = 1000000
    else:
        digits_to_round = int(5-np.log10(tick_size))
        divisor = 100000

    tick_labels = [ str(round(i/float(divisor), digits_to_round))+scale
                   for i in ticks ]

    if(loc_coordinates == "up"):
        plt.plot([start, end],
                 [0, 0],
                 linestyle="-",
                 color=color,
                 linewidth=1)
    else:
        plt.plot([start, end],
                 [0.3, 0.3],
                 linestyle="-",
                 color=color,
                 linewidth=1)


    if(revert_coordinates):
        ticks = [ start + end-i for i in ticks ]
        ticks.reverse()
        tick_labels.reverse()
    print(tick_labels)
    for i in range(len(ticks)):
        if(loc_coordinates == "up"):
            plt.plot([ticks[i], ticks[i]],
                     [0., .3],
                     linestyle="-",
                     color=color,
                     linewidth=1)
            plt.text(ticks[i],
                     .4,
                     tick_labels[i],
                     horizontalalignment="center",
                     verticalalignment="bottom",
                     fontsize=5,
                     color=color,
                     rotation=rotation)
        else:
            plt.plot([ticks[i], ticks[i]],
                     [.3, .0],
                     linestyle="-",
                     color=color,
                     linewidth=1)
            plt.text(ticks[i],
                     -.1,
                     tick_labels[i],
                     horizontalalignment="center",
                     fontsize=5,
                     color=color,
                     verticalalignment="top",
                     rotation=rotation)

    plt.xlim([start, end])
    plt.yticks([], [])
    if(loc_coordinates == "up"):
        plt.ylim([-.1, .8])
    else:
        plt.ylim([-1.5, .3])
    plt.xticks([], [])
    ax.spines["bottom"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)

