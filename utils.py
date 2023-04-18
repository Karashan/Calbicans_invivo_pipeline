# Author: Xiang Zhang
# Description: all utility functions for generating LFC values given C. albicans in vivo data and make QC plots
# support pre-processing and post-processing pipeline

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from functools import reduce
from statsmodels.stats.multitest import fdrcorrection


def create_directories(direc_list):
    for direc_str in direc_list:
        if not os.path.exists(direc_str):
            os.makedirs(direc_str)


def plot_inoc_qc(df, inoc_list, up_list, dn_list, plates, plate_col, read_thre_good=50, read_thre_bad=0,
                 inoc_prefix='SI_', output_direc='./', plt_shape=None):
    if not os.path.exists(output_direc):
        os.makedirs(output_direc)
    plot_fraction_strains_with_good_inoc(df, inoc_list=inoc_list, plates=plates, plate_col=plate_col,
                                         read_thre=read_thre_good, output_direc=output_direc + inoc_prefix, plt_shape=plt_shape)
    plot_fraction_strains_with_no_inoc(df, inoc_list=inoc_list, plates=plates, plate_col=plate_col,
                                       read_thre=read_thre_good, output_direc=output_direc + inoc_prefix, plt_shape=plt_shape)
    plot_fraction_dropout_strains(df, up_list, dn_list, plates=plates, plate_col=plate_col, read_thre=read_thre_bad,
                                  inoc_list=inoc_list, mask_val=read_thre_good,
                                  output_direc=output_direc + inoc_prefix + '_dropout_', plt_shape=plt_shape)


def df_BH_correct_pval(df, col, pval_thre=0.05):
    df['FDR'] = np.nan
    df_idx = df[df[col].notnull()]
    (reject_list, pval_bh) = fdrcorrection(df_idx[col], alpha=pval_thre, method='indep', is_sorted=False)
    df.loc[df_idx.index, 'FDR'] = pval_bh
    return df


# Read count analysis:
def find_nlargest_strains_in_plate(df_raw, plate_num, mice_list, nlargest=3, threshold=0, plate_col='plate',
                                   orf19_col='orf19', common_col='Common'):
    df_raw_plate = df_raw[df_raw[plate_col] == plate_num]
    mice_entry_list_ID = []
    mice_list_dup = []

    mice_top1_list = []
    mice_top2_list = []
    mice_top3_list = []
    orf19_list = []
    common_name_list = []
    for m in mice_list:
        df_plate_m = df_raw_plate[m]
        plate_read_sum = df_plate_m.sum()

        if plate_read_sum <= 0:
            continue
        else:
            df_nlargest_mice = df_plate_m.nlargest(nlargest)

            mice_percentage = df_nlargest_mice.values / plate_read_sum
            mice_percentage = np.around(mice_percentage * 100, 2)
            if np.sum(mice_percentage) > threshold:
                mice_entry = ''
                orf19_entry = ''
                name_entry = ''
                mice_list_dup.append(m)
                mice_top1_list.append(mice_percentage[0])
                mice_top2_list.append(np.around(mice_percentage[0] + mice_percentage[1], 2))
                mice_top3_list.append(np.around(mice_percentage.sum(), 2))
                for i in range(nlargest):
                    if i == nlargest - 1:
                        mice_entry += (df_nlargest_mice.index[i] + ' (' + str(mice_percentage[i]) + '%)')
                        orf19_entry += str(df_raw_plate.loc[df_nlargest_mice.index[i]][orf19_col])
                        name_entry += str(df_raw_plate.loc[df_nlargest_mice.index[i]][common_col])
                    else:
                        mice_entry += (df_nlargest_mice.index[i] + ' (' + str(mice_percentage[i]) + '%); ')
                        orf19_entry += str(df_raw_plate.loc[df_nlargest_mice.index[i]][orf19_col]) + '; '
                        name_entry += str(df_raw_plate.loc[df_nlargest_mice.index[i]][common_col]) + '; '

                mice_entry_list_ID.append(mice_entry)
                orf19_list.append(orf19_entry)
                common_name_list.append(name_entry)

    plate_list = [plate_num] * len(mice_list_dup)
    return pd.DataFrame(
        data={plate_col: plate_list, 'Mouse': mice_list_dup, 'Top1(%)': mice_top1_list, 'Top2(%)': mice_top2_list,
              'Top3(%)': mice_top3_list, 'Top3_strains': mice_entry_list_ID, 'orf19_in_order': orf19_list,
              'name_in_order': common_name_list}, index=None)


# Functions regarding getting LFC values
def normalize_reads(df, scaling_factor=1e6, pseudocount=1):
    return np.log2((df / df.sum()) * scaling_factor + pseudocount)


def get_masked_inoc_strains_separate(df, inoc_list, df_all, mask_val, output_direc):
    for inoc in inoc_list:
        idx_nan = df[df[inoc].isnull()].index
        df_nan = df_all.loc[idx_nan]
        df_nan.to_csv(output_direc + "/filter_" + inoc + "_less_than_" + str(mask_val) + ".csv")


def get_masked_inoc_strains_intersect(df, inoc_list, df_all, mask_val, output_direc):
    all_nan_idx = []
    for inoc in inoc_list:
        idx_nan = df[df[inoc].isnull()].index.tolist()
        if len(idx_nan) > 0:
            all_nan_idx.append(idx_nan)
    if len(all_nan_idx) > 0:
        all_nan_idx = list(reduce(set.intersection, [set(item) for item in all_nan_idx]))
    df_nan = df_all.loc[all_nan_idx]
    df_nan.to_csv(output_direc + "/filter_all_inocs_less_than_" + str(mask_val) + ".csv")


def get_corr_replicate_plots(df, output_direc, name='lfc', corr_method='pearson'):
    df_corr = df.corr(method=corr_method)
    df_corr.to_csv(output_direc + 'replicates_correlation_' + name + '.csv')
    pd.plotting.scatter_matrix(df.dropna(axis=1, how='all'), alpha=1.0, hist_kwds={'bins': 30}, figsize=(35, 25))
    plt.savefig(output_direc + 'replicates_plot_' + name + '.pdf', dpi=300)
    plt.close()


def plot_fraction_strains_with_good_inoc(df, inoc_list, plates=None, plate_col='plate', read_thre=50,
                                         output_direc='fraction_plot/', plt_shape=None):
    if plates is None:
        plates = df[plate_col].unique()
    if len(inoc_list) != 2:
        print(
            "To plot fraction of strains with qualified inoculum reads, please give an inoc_list with length of 2 (usually 1 UP & 1 DN tag).")
        return None

    plate_list = []
    for plate in plates:
        df_plate = df[df[plate_col] == plate]
        plate_inoc_list = []
        for inoc in inoc_list:
            df_plate_inoc = df_plate[inoc]
            df_plate_inoc_thre = df_plate_inoc[df_plate_inoc > read_thre]
            strain_fraction = float(df_plate_inoc_thre.shape[0] / df_plate.shape[0])
            plate_inoc_list.append(round(strain_fraction, 4))
        plate_list.append(plate_inoc_list)

    plate_array = np.array(plate_list).T * 100

    # Grouped bar chart
    # labels = plates
    x = np.arange(len(plates))
    width = 0.35  # the width of the bars

    if plt_shape != None:
        fig, ax = plt.subplots(figsize=plt_shape)
    else:
        fig, ax = plt.subplots()
    ax.set_ylim(top=105)
    rects1 = ax.bar(x - width / 2, plate_array[0], width, label=inoc_list[0])
    rects2 = ax.bar(x + width / 2, plate_array[1], width, label=inoc_list[1])

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Percentage of strains with more than ' + str(read_thre) + " reads in inoculum (%)")
    # ax.set_title('')
    ax.set_xticks(x)
    ax.set_xticklabels(plates)
    ax.legend()

    ax.bar_label(rects1, padding=3)
    ax.bar_label(rects2, padding=3)

    fig.tight_layout()

    plt.savefig(output_direc + "_fraction_strains_inoc_more_than_" + str(read_thre) + ".pdf", dpi=300)
    plt.close()


def plot_fraction_strains_with_no_inoc(df, inoc_list, plates=None, plate_col='plate', read_thre=50,
                                       output_direc='fraction_plot/', plt_shape=None):
    if plates is None:
        plates = df[plate_col].unique()
    if len(inoc_list) != 2:
        print("To plot fraction of strains with qualified inoculum reads, please give an inoc_list with length of 2 (usually 1 UP & 1 DN tag).")
        return None

    plate_list = []
    for plate in plates:
        df_plate = df[df[plate_col] == plate]
        plate_inoc_list = []
        # for inoc in inoc_list:
        df_plate_inoc_up = df_plate[inoc_list[0]]
        df_plate_inoc_dn = df_plate[inoc_list[1]]
        df_plate_inoc_up_thre = df_plate_inoc_up[df_plate_inoc_up < read_thre]
        df_plate_inoc_dn_thre = df_plate_inoc_dn[df_plate_inoc_dn < read_thre]
        no_coverage_idx = np.intersect1d(df_plate_inoc_up_thre.index, df_plate_inoc_dn_thre.index)
        strain_fraction = float(len(no_coverage_idx) / df_plate.shape[0])
        plate_inoc_list.append(round(strain_fraction, 4))
        plate_list.append(plate_inoc_list)

    plate_array = np.array(plate_list).T * 100

    # Grouped bar chart
    # labels = plates
    x = np.arange(len(plates))
    width = 0.35  # the width of the bars

    if plt_shape != None:
        fig, ax = plt.subplots(figsize=plt_shape)
    else:
        fig, ax = plt.subplots()
    ax.set_ylim(top=105)
    rects1 = ax.bar(x, plate_array[0], width, label=inoc_list[0][:-3])
    # rects2 = ax.bar(x + width/2, plate_array[1], width, label=inoc_list[1])

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Percentage of strains with less than ' + str(read_thre) + " reads in inoculum (%)")
    # ax.set_title('')
    ax.set_xticks(x)
    ax.set_xticklabels(plates)
    ax.legend()

    ax.bar_label(rects1, padding=3)
    # ax.bar_label(rects2, padding=3)

    fig.tight_layout()

    plt.savefig(output_direc + "_fraction_strains_inoc_less_than_" + str(read_thre) + ".pdf", dpi=300)
    plt.close()


#  Old version where we do not exclude any background strains
def plot_fraction_dropout_strains_old(df, rep_up_list, rep_dn_list, plates=None, plate_col='plate', read_thre=0,
                                  output_direc='fraction_plot/', plt_shape=None):
    if plates is None:
        plates = df[plate_col].unique()

    plate_list = []
    for plate in plates:
        df_plate = df[df[plate_col] == plate]
        plate_read_list = []
        # print(plate)
        for i in range(len(rep_up_list)):  # assume same length with dn list
            df_plate_up = df_plate[rep_up_list[i]]
            df_plate_dn = df_plate[rep_dn_list[i]]
            df_plate_up_thre = df_plate_up[df_plate_up <= read_thre]
            df_plate_dn_thre = df_plate_dn[df_plate_dn <= read_thre]
            no_coverage_idx = np.intersect1d(df_plate_up_thre.index, df_plate_dn_thre.index)
            strain_fraction = float(len(no_coverage_idx) / df_plate.shape[0])
            plate_read_list.append(round(strain_fraction, 4) * 100)
        # print(plate_read_list)
        plate_list.append(plate_read_list)
    col_names = [rep.replace('_UP', '') for rep in rep_up_list]
    df_result = pd.DataFrame(data=plate_list, index=plates, columns=col_names)

    fig, ax = plt.subplots()
    ax = df_result.plot(y=col_names, kind="bar", rot=0, figsize=plt_shape)
    ax.set_ylim(top=105)
    for i in range(len(col_names)):
        ax.bar_label(ax.containers[i], padding=3, size=3)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel("Percentage of dropout strains (0 read in UP&DN) (%)")
    fig.tight_layout()

    plt.savefig(output_direc + "fraction_dropout_strains_in_reps.pdf", dpi=300)
    plt.close()


# New version where we calculate dropout using the background excluding the strains that could not pass QC
def plot_fraction_dropout_strains(df, rep_up_list, rep_dn_list, plates=None, plate_col='plate', read_thre=0,
                                  inoc_list=None, mask_val=50, output_direc='fraction_plot/', plt_shape=None):

    if plates is None:
        plates = df[plate_col].unique()

    all_nan_idx = []
    if inoc_list is None:
        print("To plot the fraction of dropout strains with correct background, please give an inoc_list with length of 2 (usually 1 UP & 1 DN tag).")
    else:
        for inoc in inoc_list:
            idx_nan = df[df[inoc] < mask_val].index.tolist()
            if len(idx_nan) > 0:
                all_nan_idx.append(idx_nan)
        if len(all_nan_idx) > 0:
            all_nan_idx = list(reduce(set.intersection, [set(item) for item in all_nan_idx]))
        print("Exclude", len(all_nan_idx), "strains that could not pass the", str(mask_val),
              "reads threshold for either inoculum in all plates.\n")
    df_ = df.drop(all_nan_idx, axis=0)

    plate_list = []
    for plate in plates:
        df_plate = df_[df_[plate_col] == plate]
        plate_read_list = []
        # print(plate)
        for i in range(len(rep_up_list)):  # assume same length with dn list
            df_plate_up = df_plate[rep_up_list[i]]
            df_plate_dn = df_plate[rep_dn_list[i]]
            df_plate_up_thre = df_plate_up[df_plate_up <= read_thre]
            df_plate_dn_thre = df_plate_dn[df_plate_dn <= read_thre]
            no_coverage_idx = np.intersect1d(df_plate_up_thre.index, df_plate_dn_thre.index)
            strain_fraction = float(len(no_coverage_idx) / df_plate.shape[0])
            plate_read_list.append(round(strain_fraction, 4) * 100)
        # print(plate_read_list)
        plate_list.append(plate_read_list)
    col_names = [rep.replace('UP', '') for rep in rep_up_list]
    col_names = [col.strip() for col in col_names]
    df_result = pd.DataFrame(data=plate_list, index=plates, columns=col_names)

    fig, ax = plt.subplots()
    ax = df_result.plot(y=col_names, kind="bar", rot=0, figsize=plt_shape)
    ax.set_ylim(top=105)
    for i in range(len(col_names)):
        ax.bar_label(ax.containers[i], padding=3, size=3)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel("Percentage of dropout strains (0 read in UP&DN) (%)")
    fig.tight_layout()

    plt.savefig(output_direc + "fraction_dropout_strains_in_reps.pdf", dpi=300)
    plt.close()


def get_LFC_with_mask_singleInoc(df_all, name_col_list, inoc_up, inoc_dn, mouse_up_list, mouse_dn_list,
                                 plate_col='plate', plate=None, output_lfc_direc=None, mask_val=50, make_QC=True):
    if plate != None:
        df_all = df_all[df_all[plate_col] == plate]

    # Data preparation
    inoc_list = [inoc_up, inoc_dn]
    mouse_all_list = mouse_up_list + mouse_dn_list
    mouse_all_lfc_list = [m + '_LFC' for m in mouse_all_list]

    # Filter: replace any inoc reads < mask_val (e.g. 50) with NaN
    df_inoc = df_all[inoc_list]
    df_inoc_masked = df_inoc.mask(df_inoc < mask_val)
    df_inoc_masked[mouse_all_list] = df_all[mouse_all_list]

    # Normalize and get Log2 values
    df_inoc_masked_normalized = normalize_reads(df_inoc_masked)

    # Get LFC
    for m_up in mouse_up_list:
        m_up_lfc = m_up + '_LFC'
        df_inoc_masked_normalized[m_up_lfc] = df_inoc_masked_normalized[m_up] - df_inoc_masked_normalized[inoc_up]
    for m_dn in mouse_dn_list:
        m_dn_lfc = m_dn + '_LFC'
        df_inoc_masked_normalized[m_dn_lfc] = df_inoc_masked_normalized[m_dn] - df_inoc_masked_normalized[inoc_dn]

    for name_col in name_col_list:
        df_inoc_masked_normalized[name_col] = df_all[name_col]

    df_inoc_masked_normalized = df_inoc_masked_normalized[
        name_col_list + inoc_list + mouse_all_list + mouse_all_lfc_list]

    if output_lfc_direc != None:
        df_inoc_masked_normalized[name_col_list + mouse_all_lfc_list].to_csv(output_lfc_direc + '_LFC.csv')

    # If complete data was given with multiple plates, we iterate through plates
    if make_QC:
        qc_direc = output_lfc_direc + '_QC/'
        if not os.path.exists(qc_direc):
            os.makedirs(qc_direc)

        # Get the strains that are filtered due to inoc reads < mask_val (e.g. 50)
        get_masked_inoc_strains_intersect(df_inoc_masked_normalized, inoc_list, df_all, mask_val, qc_direc)

        # Get the correlation table for both the raw reads and the LFC values
        # get_corr_replicate_plots(df_all[mouse_all_list], qc_direc, name='raw_reads_all')
        # get_corr_replicate_plots(df_inoc_masked_normalized[mouse_all_lfc_list], qc_direc, name='lfc_all')

        # Separate QC plots for each plate individually if there are multiple plates in the input df
        plates = df_all[plate_col].unique().tolist()
        for p in plates:
            df_p = df_all[df_all[plate_col] == p]
            p_direc = qc_direc + p + '/'
            if not os.path.exists(p_direc):
                os.makedirs(p_direc)
            get_corr_replicate_plots(df_p[mouse_all_list], p_direc, name='raw_reads_' + p)
            get_corr_replicate_plots(df_inoc_masked_normalized.loc[df_p.index, mouse_all_lfc_list], p_direc,
                                     name='lfc_' + p)

    return df_inoc_masked_normalized


# Assume same tags (UP/DN) are input
def get_chosen_LFC_from_2inoc(inoc1, inoc2, df, mouse):
    lfc_col = mouse + '_LFC'

    df_inoc_m = df[[inoc1, inoc2, mouse]]
    df_inoc_m['inoc1_LFC'] = df_inoc_m[mouse] - df_inoc_m[inoc1]
    df_inoc_m['inoc2_LFC'] = df_inoc_m[mouse] - df_inoc_m[inoc2]
    inoc1_notnull_idx = df_inoc_m['inoc1_LFC'][df_inoc_m['inoc1_LFC'].notnull()].index
    inoc2_notnull_idx = df_inoc_m['inoc2_LFC'][df_inoc_m['inoc2_LFC'].notnull()].index

    df_inoc_m[lfc_col] = np.nan
    # Cover all possible LFCs with both inocs
    df_inoc_m.loc[inoc1_notnull_idx, lfc_col] = df_inoc_m.loc[inoc1_notnull_idx, 'inoc1_LFC']
    df_inoc_m.loc[inoc2_notnull_idx, lfc_col] = df_inoc_m.loc[inoc2_notnull_idx, 'inoc2_LFC']
    # Change the overlap part to the average LFC of both
    overlap_idx = np.intersect1d(inoc1_notnull_idx, inoc2_notnull_idx)
    df_inoc_m.loc[overlap_idx, lfc_col] = (df_inoc_m.loc[overlap_idx, 'inoc1_LFC'] + df_inoc_m.loc[
        overlap_idx, 'inoc2_LFC']) / 2

    return df_inoc_m[lfc_col]


def get_LFC_with_mask_twoInoc(df_all, name_col_list, inoc_up_list, inoc_dn_list, mouse_up_list, mouse_dn_list,
                              plate_col='plate', plate=None, output_lfc_direc=None, mask_val=50, make_QC=True):
    # Focus on the only re-arrayed plate for accurate normalization
    # Ideally, this would be done across each sequencing lane
    if plate != None:
        df_all = df_all[df_all[plate_col] == plate]

    # Data preparation
    # The inoc_up_list and inoc_dn_list must each contain two UP or DN tags from different inoculums
    inoc_list = inoc_up_list + inoc_dn_list
    mouse_all_list = mouse_up_list + mouse_dn_list
    mouse_all_lfc_list = [m + '_LFC' for m in mouse_all_list]

    # Filter: replace any inoc reads < mask_val (e.g. 50) with NaN
    df_inoc = df_all[inoc_list]
    df_inoc_masked = df_inoc.mask(df_inoc < mask_val)
    df_inoc_masked[mouse_all_list] = df_all[mouse_all_list]

    # Normalize and get Log2 values
    df_inoc_masked_normalized = normalize_reads(df_inoc_masked)

    # Get LFC
    for m_up in mouse_up_list:
        m_up_lfc = m_up + '_LFC'
        df_inoc_masked_normalized[m_up_lfc] = get_chosen_LFC_from_2inoc(inoc_up_list[0], inoc_up_list[1],
                                                                        df_inoc_masked_normalized, m_up)
    for m_dn in mouse_dn_list:
        m_dn_lfc = m_dn + '_LFC'
        df_inoc_masked_normalized[m_dn_lfc] = get_chosen_LFC_from_2inoc(inoc_dn_list[0], inoc_dn_list[1],
                                                                        df_inoc_masked_normalized, m_dn)

    for name_col in name_col_list:
        df_inoc_masked_normalized[name_col] = df_all[name_col]

    df_inoc_masked_normalized = df_inoc_masked_normalized[
        name_col_list + inoc_list + mouse_all_list + mouse_all_lfc_list]

    # Create QC folder regardless of the setting of make_QC
    qc_direc = output_lfc_direc + '_QC'
    if not os.path.exists(qc_direc):
        os.makedirs(qc_direc, exist_ok=True)
    if output_lfc_direc != None:
        df_inoc_masked_normalized[name_col_list + mouse_all_lfc_list].to_csv(output_lfc_direc + '.csv')

    if make_QC:
        # Get the strains that are filtered due to inoc reads < mask_val (e.g. 50) across all inoculums
        get_masked_inoc_strains_intersect(df_inoc_masked_normalized, inoc_list, df_all, mask_val, qc_direc)

        # If complete data was given with multiple plates, here we should iterate through plates 
        for p in df_all[plate_col].unique():
            plate_direc = qc_direc + '/' + p.replace(' ', '')
            if not os.path.exists(plate_direc):
                os.makedirs(plate_direc)
            plate_idx = df_all[df_all[plate_col] == p].index
            # Get the correlation table for both the raw reads and the LFC values
            get_corr_replicate_plots(df_all.loc[plate_idx, inoc_list + mouse_all_list], plate_direc, name='raw_reads')
            get_corr_replicate_plots(df_inoc_masked_normalized.loc[plate_idx, inoc_list + mouse_all_list], plate_direc,
                                     name='log2')
            get_corr_replicate_plots(df_inoc_masked_normalized.loc[plate_idx, mouse_all_lfc_list], plate_direc,
                                     name='lfc')

    return df_inoc_masked_normalized


def get_topN_dom_strains(df_rearray_feature, plate_list, mice_list, output_direc, topN=3, thre=0, plate_col='Plate',
                         orf19_col='orf19', common_col='Common'):
    list_df_plate_strain_stat = []
    for plate_num in plate_list:
        list_df_plate_strain_stat.append(find_nlargest_strains_in_plate(df_rearray_feature, plate_num,
                                                                        mice_list=mice_list, nlargest=topN,
                                                                        threshold=thre,
                                                                        plate_col=plate_col, orf19_col=orf19_col,
                                                                        common_col=common_col))
    df_plate_strain_stat_concat = pd.concat(list_df_plate_strain_stat)
    if output_direc != None:
        df_plate_strain_stat_concat.to_csv(output_direc, index=None)
    return df_plate_strain_stat_concat


# Correlation analysis
def produce_corr_table_across_mice(df_lfc, input_plate_list, sub_lfc_list, df_dominate, tail='_LFC', method_='pearson',
                                   plate_col='plate'):
    pcc_list = []
    plate_list = []
    plate_M_list = []
    for plate in input_plate_list:
        df_lfc_plate_corr_pcc = df_lfc.loc[df_lfc[df_lfc[plate_col] == plate].index, sub_lfc_list].corr(method=method_)
        df_dominate_plate = df_dominate[df_dominate[plate_col] == plate]

        for row in df_dominate_plate.iterrows():
            m = row[1]['Mouse']
            m_list = df_dominate_plate['Mouse'].tolist()
            m_list.remove(m)
            pcc_sub_list = []
            for m_other in m_list:
                pcc_sub = df_lfc_plate_corr_pcc.loc[m + tail, m_other + tail]
                pcc_sub_list.append(pcc_sub)
            plate_list.append(plate)
            plate_M_list.append(m + tail)
            pcc_list.append(np.mean(pcc_sub_list))
    df = pd.DataFrame({plate_col: plate_list, 'Mouse': plate_M_list, 'mean_pcc_with_others': pcc_list})
    return df


def plot_avgCorr_vs_pcc(df_top3, df_lfc, plate_list, filename, plate_col='plate'):
    for plate in plate_list:
        plt.scatter(df_top3[df_top3[plate_col] == plate]['Top3(%)'],
                    df_lfc[df_lfc[plate_col] == plate]['mean_pcc_with_others'])
    plt.legend(plate_list, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('Sum of Top3 strains (%)')
    plt.ylabel('Average correlation between the LFC of mouse and other mice')
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()


def get_avgCorr_vs_sumTop3(lfc_direc, dom_up_direc, dom_dn_direc, m_lfc_list_up, m_lfc_list_dn, output_direc,
                           tail_='_LFC', corr_method='pearson', plate_col='plate', sort_cols=['plate', 'Mouse']):
    # 2. Get LFC values
    df_LFC = pd.read_csv(lfc_direc, index_col=0)
    plate_list = df_LFC[plate_col].unique().tolist()

    # 3. Get dominating Top 3 information
    df_dom_UP = pd.read_csv(dom_up_direc, index_col=None)
    df_dom_DN = pd.read_csv(dom_dn_direc, index_col=None)

    # 4. Get the correlation table
    df_UP_lfc = produce_corr_table_across_mice(df_LFC, plate_list, m_lfc_list_up, df_dom_UP, tail=tail_,
                                               method_=corr_method, plate_col=plate_col)
    df_UP_lfc = df_UP_lfc.sort_values(sort_cols)
    df_DN_lfc = produce_corr_table_across_mice(df_LFC, plate_list, m_lfc_list_dn, df_dom_DN, tail=tail_,
                                               method_=corr_method, plate_col=plate_col)
    df_DN_lfc = df_DN_lfc.sort_values(sort_cols)

    df_dom_UP = df_dom_UP.sort_values(sort_cols)
    df_dom_DN = df_dom_DN.sort_values(sort_cols)

    df_UP_lfc.set_index(plate_col).to_csv(output_direc + '_UP.txt', sep='\t')
    df_DN_lfc.set_index(plate_col).to_csv(output_direc + '_DN.txt', sep='\t')

    # PCC between LFC
    plot_avgCorr_vs_pcc(df_dom_UP, df_UP_lfc, plate_list, filename=output_direc + '_UP.pdf', plate_col=plate_col)
    plot_avgCorr_vs_pcc(df_dom_DN, df_DN_lfc, plate_list, filename=output_direc + '_DN.pdf', plate_col=plate_col)


# Variance normalization
def normalize_variance(df_LFC, name_LFC_list, plate_col='plate'):
    plate_list = df_LFC[plate_col].unique()
    df_p_list = []
    for p in plate_list:
        df_LFC_p = df_LFC[df_LFC[plate_col] == p]
        # Goal: let every set of experiments be normalized to get a consistent variance
        df_p_t_list = []
        for tag in name_LFC_list:
            df_LFC_p_t = df_LFC_p[tag]
            alpha_R_LFC = df_LFC_p_t.mean(axis=1, skipna=True).std(skipna=True)
            scalars_ar_am = alpha_R_LFC / df_LFC_p_t.std(skipna=True).values  # a_r/a_m
            df_LFC_p_t_norm = df_LFC_p_t * scalars_ar_am
            df_p_t_list.append(df_LFC_p_t_norm)
        # Concatanate the list for this plate
        df_p = pd.concat(df_p_t_list, axis=1)
        plate_sub_list = [p] * len(df_LFC_p)
        df_p[plate_col] = plate_sub_list
        df_p_list.append(df_p)
    df_LFC_norm = pd.concat(df_p_list)
    return df_LFC_norm


def generate_modt_matrices(df_LFC_rm, column_names, UPorDN, output_direc, plate_col='plate', file_head='SI_LFC_plate'):
    if not os.path.exists(output_direc):
        os.makedirs(output_direc)
    for plate in df_LFC_rm[plate_col].unique():
        plate_orf = df_LFC_rm[df_LFC_rm[plate_col] == plate].index
        df_plate = df_LFC_rm.loc[plate_orf, column_names]
        df_plate = df_plate.dropna(axis=1, how='all')
        df_plate = df_plate.dropna(axis=0, how='all')
        plate_orf_new = df_plate.index
        vec_collect = []
        for row in df_plate.iterrows():
            idx = row[0]
            vec_first = row[1].values
            vec_others = df_plate.drop(idx).stack().values
            vec_temp = np.append(vec_first, vec_others)
            vec_collect.append(vec_temp)
        df_plate_matrix = pd.DataFrame(np.array(vec_collect))
        df_plate_matrix.index = plate_orf_new
        df_plate_matrix.to_csv(
            (output_direc + file_head + str(plate[-2:]) + '_' + UPorDN + '_' + str(df_plate.shape[1]) + 'rep' + '.txt'),
            sep='\t', header=None)
    return df_plate_matrix


def generate_merged_modt_matrices(df_LFC_rm_norm_direc, m_lfc_list_up, m_lfc_list_dn, output_direc, plate_col='plate',
                                  file_head='SI_LFC_plate', minimal_cutoff=3):
    if not os.path.exists(output_direc):
        os.makedirs(output_direc)
    df_LFC_rm = pd.read_csv(df_LFC_rm_norm_direc, sep='\t', index_col=0)
    column_names = []
    # Assume UP and DN LFC lists are always parallel with the same length
    for i in range(len(m_lfc_list_up)):
        column_names.append(m_lfc_list_up[i])
        column_names.append(m_lfc_list_dn[i])

    for plate in df_LFC_rm[plate_col].unique():
        plate_orf = df_LFC_rm[df_LFC_rm[plate_col] == plate].index
        df_plate = df_LFC_rm.loc[plate_orf, column_names]
        # df_plate = df_plate.dropna(axis=1, how='all')
        # Remove rows with all NaNs, but leave all columns there for structual need
        df_plate = df_plate.dropna(axis=0, how='all')
        # plate_orf_new = df_plate.index
        vec_collect = []
        idx_collect = []
        for row in df_plate.iterrows():
            idx = row[0]
            vec_first = row[1].values
            if np.count_nonzero(~np.isnan(vec_first)) < minimal_cutoff:
                continue
            vec_others = df_plate.drop(idx).stack(dropna=False).values
            vec_temp = np.append(vec_first, vec_others)
            vec_collect.append(vec_temp)
            idx_collect.append(idx)
        df_plate_matrix = pd.DataFrame(np.array(vec_collect))
        if df_plate_matrix.shape[0] == 0:
            continue
        df_plate_matrix.index = idx_collect
        df_plate_matrix.to_csv((output_direc + file_head + str(plate.replace(' ', '')) + '_' + str(
            int(df_plate.shape[1] / 2)) + 'rep' + '.txt'), sep='\t', header=None)
    return df_plate_matrix


# Update: changed filter criteria to include both tags (mouse-based)
def get_var_norm_mod_t_input_with_filter(dataset_type, lfc_direc, dom_direc, avg_up_direc, avg_dn_direc,
                                         m_lfc_list_up, m_lfc_list_dn, lfc_rm_norm_output_direc,
                                         lfc_rm_norm_mod_output_direc, tail_='_LFC', plate_col='plate', dom_thre=90,
                                         corr_thre=0.5, up_tag='UP', dn_tag='DN'):
    # 1. Get LFC values
    sep_ = ','
    if '.txt' in lfc_direc:
        sep_ = '\t'
    df_LFC = pd.read_csv(lfc_direc, sep=sep_, index_col=0)
    if plate_col not in df_LFC.columns:
        print("Warning: LFC file does not have plate info!")

    print("Dominating strains threshold (%):", dom_thre)
    print("Average correlation threshold:", corr_thre)

    # 2. Get dominating Top 3 information, remove the top3 > threshold
    df_dom = pd.read_csv(dom_direc, index_col=0)
    df_dom_top3thre = df_dom[df_dom['Top3(%)'] > dom_thre]
    for row in df_dom_top3thre.iterrows():
        plate_ = row[0]
        mouse_ = row[1]['Mouse']
        mouse_ = mouse_ + tail_
        if mouse_ not in df_LFC.columns:  # double safe check
            print(mouse_ + " not in this plate")
            continue

        plate_idx = df_LFC[df_LFC[plate_col] == plate_].index
        if len(plate_idx) > 0:
            df_LFC.loc[plate_idx, mouse_] = np.nan
            print("Top3: Excluding " + plate_ + ' ' + mouse_)
            # Deal with the other tag
            if up_tag in mouse_:
                mouse_other_tag = mouse_.replace(up_tag, dn_tag)
            elif dn_tag in mouse_:
                mouse_other_tag = mouse_.replace(dn_tag, up_tag)
            else:
                print("Error when getting the other tag of the mouse")
            df_LFC.loc[plate_idx, mouse_other_tag] = np.nan
            print("Top3: Excluding " + plate_ + ' ' + mouse_other_tag)

    # 3. Remove flagged mice & plate from the correlation analysis
    # SI flag list (threshold: 0.1 correlation coefficient):
    # UP
    df_avgCorr_up = pd.read_csv(avg_up_direc, sep='\t', index_col=None)
    # DN
    df_avgCorr_dn = pd.read_csv(avg_dn_direc, sep='\t', index_col=None)

    df_avgCorr_up_thre = df_avgCorr_up[df_avgCorr_up['mean_pcc_with_others'] < corr_thre]
    df_avgCorr_dn_thre = df_avgCorr_dn[df_avgCorr_dn['mean_pcc_with_others'] < corr_thre]

    for row in df_avgCorr_up_thre.iterrows():
        plate_ = row[1][plate_col]
        mouse_ = row[1]['Mouse']
        print("Corr: Excluding " + plate_ + ' ' + mouse_)
        df_LFC.loc[df_LFC[df_LFC[plate_col] == plate_].index, mouse_] = np.nan
        # Deal with the other tag
        if up_tag in mouse_:
            mouse_other_tag = mouse_.replace(up_tag, dn_tag)
        elif dn_tag in mouse_:
            mouse_other_tag = mouse_.replace(dn_tag, up_tag)
        else:
            print("Error when getting the other tag of the mouse")
        print("Corr: Excluding " + plate_ + ' ' + mouse_other_tag)
        df_LFC.loc[df_LFC[df_LFC[plate_col] == plate_].index, mouse_other_tag] = np.nan

    for row in df_avgCorr_dn_thre.iterrows():
        plate_ = row[1][plate_col]
        mouse_ = row[1]['Mouse']
        print("Corr: Excluding " + plate_ + ' ' + mouse_)
        df_LFC.loc[df_LFC[df_LFC[plate_col] == plate_].index, mouse_] = np.nan
        # Deal with the other tag
        if up_tag in mouse_:
            mouse_other_tag = mouse_.replace(up_tag, dn_tag)
        elif dn_tag in mouse_:
            mouse_other_tag = mouse_.replace(dn_tag, up_tag)
        else:
            print("Error when getting the other tag of the mouse")
        print("Corr: Excluding " + plate_ + ' ' + mouse_other_tag)
        df_LFC.loc[df_LFC[df_LFC[plate_col] == plate_].index, mouse_other_tag] = np.nan

    # Save the filtered but unnormalized LFC values
    df_LFC.to_csv(lfc_rm_norm_output_direc.replace('_var_norm', ''), sep='\t')

    # 4. Output a file that includes the normalization(Calculate M_hat) and before normalization
    df_LFC_rm_norm = normalize_variance(df_LFC, [m_lfc_list_up, m_lfc_list_dn], plate_col=plate_col)
    df_LFC_rm_norm.to_csv(lfc_rm_norm_output_direc, sep='\t')

    # 5. Generate the input matrices for moderated t-test
    # df_up_var = generate_modt_matrices(df_LFC_rm_norm, m_lfc_list_up, 'UP',
    #                                    output_direc=lfc_rm_norm_mod_output_direc + '_UP/', file_head=dataset_type,
    #                                    plate_col=plate_col)
    # df_dn_var = generate_modt_matrices(df_LFC_rm_norm, m_lfc_list_dn, 'DN',
    #                                    output_direc=lfc_rm_norm_mod_output_direc + '_DN/', file_head=dataset_type,
    #                                    plate_col=plate_col)
    generate_merged_modt_matrices(lfc_rm_norm_output_direc, m_lfc_list_up, m_lfc_list_dn,
                                  output_direc=lfc_rm_norm_mod_output_direc + '/',
                                  plate_col=plate_col, file_head=dataset_type)


def split_str_list(x, spl='|'):
    if x is not np.nan:
        x = x.split(spl)
    else:
        x = np.nan
    return x


# Can add relevant gene info including ASSEMBLY22_ID, common name, description, and S. cer ortholog to a dataframe with C.albicans orf19 as index
# Make sure downloading the latest files from CGD
def add_common_name(df,
                    df_maporf19_direc='/Users/zhangxiang/Documents/Research_CB/GeneEssentiality_CandidaAlbicans/barcode_seq/map_gene/ORF19_Assembly22_mapping.tab.txt',
                    df_descri_direc='/Users/zhangxiang/Documents/Research_CB/GeneEssentiality_CandidaAlbicans/barcode_seq/map_gene/C_albicans_SC5314_A22_current_chromosomal_feature.tab.txt',
                    df_Scer_ortho_direc='/Users/zhangxiang/Documents/Research_CB/GeneEssentiality_CandidaAlbicans/barcode_seq/map_gene/C_albicans_SC5314_S_cerevisiae_orthologs.txt'):
    input_df = df.copy(deep=True)
    df_maporf19 = pd.read_csv(df_maporf19_direc, sep='\t', index_col=0)
    df_descri = pd.read_csv(df_descri_direc, sep='\t', skiprows=8, index_col=0, header=None)
    df_Scer_ortho = pd.read_csv(df_Scer_ortho_direc, sep='\t', skiprows=8, index_col=0, header=None)

    orig_cols = input_df.columns.tolist()

    # ID and Gene name
    input_df['ASSEMBLY22_ID'] = np.nan
    input_df['GENE_NAME'] = np.nan
    overlap_orf = np.intersect1d(df_maporf19.index, input_df.index)
    input_df.loc[overlap_orf, 'ASSEMBLY22_ID'] = df_maporf19.loc[overlap_orf, 'ASSEMBLY22_ID']
    input_df.loc[overlap_orf, 'GENE_NAME'] = df_maporf19.loc[overlap_orf, 'GENE_NAME']

    # Deal with potential null values for ASSEMBLY22_ID
    df_nullID = input_df[input_df['ASSEMBLY22_ID'].isnull()]
    df_allorf = df_descri[2].apply(split_str_list)
    for orf19 in df_nullID.index:
        for i in range(len(df_allorf)):
            if df_allorf[i] is not np.nan:
                if orf19 in df_allorf[i]:
                    input_df.loc[orf19, 'ASSEMBLY22_ID'] = df_allorf.index[i]
                    print("This strain may have an updated orf19 ID:", orf19, input_df.loc[orf19, 'ASSEMBLY22_ID'])

    # Sanity check: whether clear up all NaNs or not
    null_idx = []
    df_nullID = input_df[input_df['ASSEMBLY22_ID'].isnull()]
    if df_nullID.shape[0] != 0:
        null_idx = df_nullID.index.tolist()
        print("unidentified orfs:", null_idx)

    # Description
    input_df = input_df.reset_index().set_index('ASSEMBLY22_ID')
    nonNan_idx_list = [ele for ele in input_df.index.tolist() if ele not in null_idx]
    input_df['Description'] = np.nan
    overlap_ID = np.intersect1d(df_descri.index, nonNan_idx_list)
    input_df.loc[overlap_ID, 'Description'] = df_descri.loc[overlap_ID, 10]

    # S.cer ortholog
    input_df['S.cer_Ortholog'] = np.nan
    overlap_ID = np.intersect1d(df_Scer_ortho.index, nonNan_idx_list)
    input_df.loc[overlap_ID, 'S.cer_Ortholog'] = df_Scer_ortho.loc[overlap_ID, 4]

    input_df = input_df.reset_index().set_index('orf19')
    new_columns = ['ASSEMBLY22_ID', 'GENE_NAME', 'Description', 'S.cer_Ortholog'] + orig_cols
    input_df = input_df[new_columns]
    return input_df


'''
gene names and unique identifiers for all tested mutants
descriptions of gene function from CGD
competitive index (log2[recovered/inoculum]) in each tested animal
mean competitive index across all tested animals
p-value of difference of CI from expected
Is the mutant a "hit"? (Advantage / Deffect)
'''
# Do once per dataset type (e.g. SI_GRACE)
def generate_mod_t_integrated_table(lfc_direc, m_lfc_list, mod_output_direc, output_direc, mod_effect_col='logFC',
                                    mod_meanLFC_col='mean_replicate_LFC', mod_pval_col='P.Value',
                                    hits_col='Hits(FDR<0.3,|dLFC|>2)', hits_fdr=0.3, hits_dLFC=2, plate_col='Plate',
                                    df_maporf19_direc='/Users/zhangxiang/Documents/Research_CB/GeneEssentiality_CandidaAlbicans/barcode_seq/map_gene/ORF19_Assembly22_mapping.tab.txt',
                                    df_descri_direc='/Users/zhangxiang/Documents/Research_CB/GeneEssentiality_CandidaAlbicans/barcode_seq/map_gene/C_albicans_SC5314_A22_current_chromosomal_feature.tab.txt',
                                    df_Scer_ortho_direc='/Users/zhangxiang/Documents/Research_CB/GeneEssentiality_CandidaAlbicans/barcode_seq/map_gene/C_albicans_SC5314_S_cerevisiae_orthologs.txt'
                                    ):
    df_merge = pd.read_csv(lfc_direc, sep='\t', index_col=0)  # LFC has all plate info, no need to concat
    files = os.listdir(mod_output_direc)
    df_merge['differential_LFC'] = np.nan
    df_merge['p_value'] = np.nan
    for file in files:
        if not os.path.isdir(mod_output_direc + file) and file[0] != '.':
            print(file)
            df_mod = pd.read_csv(mod_output_direc + file, index_col=0)
            if df_mod.shape[0] == 0:
                continue
            df_mod = df_mod[[mod_effect_col, mod_pval_col]]
            df_mod.columns = ['differential_LFC', 'p_value']
            idx = np.intersect1d(df_merge.index, df_mod.index)
            df_merge.loc[idx, 'differential_LFC'] = df_mod.loc[idx, 'differential_LFC']
            df_merge.loc[idx, 'p_value'] = df_mod.loc[idx, 'p_value']
    df_merge[mod_meanLFC_col] = df_merge[m_lfc_list].mean(axis=1)
    df_merge = df_BH_correct_pval(df_merge, 'p_value')
    df_merge[hits_col] = np.nan
    idx_neg_hits = df_merge[(df_merge['FDR'] < hits_fdr) & (df_merge['differential_LFC'] < -hits_dLFC)].index
    idx_pos_hits = df_merge[(df_merge['FDR'] < hits_fdr) & (df_merge['differential_LFC'] > hits_dLFC)].index
    df_merge.loc[idx_neg_hits, hits_col] = 'Defect'
    df_merge.loc[idx_pos_hits, hits_col] = 'Advantage'
    df_merge = df_merge[[plate_col, 'mean_replicate_LFC', 'differential_LFC', 'p_value', 'FDR', hits_col] + m_lfc_list]
    df_merge = add_common_name(df_merge,
                               df_maporf19_direc=df_maporf19_direc,
                               df_descri_direc=df_descri_direc,
                               df_Scer_ortho_direc=df_Scer_ortho_direc)
    df_merge.to_csv(output_direc, sep='\t')
    return df_merge
