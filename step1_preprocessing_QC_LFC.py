# Author: Xiang Zhang
# Contact: zhan6668@umn.edu
# Updated on 03/06/2023

# Step 1: With the raw reads and metatable as input,
# we generate the filtered & normalized LFC files and all the QC files for each group.
# The user needs to manually set the directories for working (where input read file is located at) and metatable.


from utils import *

if __name__ == '__main__':
    # Set commonly used parameters
    read_cutoff = 50
    dropout_cutoff = 0
    qc_plot_shape = (20, 8)

    # Set working directory
    work_dir = '/Users/zhangxiang/Documents/Research_CB/GeneEssentiality_CandidaAlbicans/barcode_seq/sequencingforseptembernoblesamples/'
    # work_dir = '/Users/zhangxiang/Documents/Research_CB/GeneEssentiality_CandidaAlbicans/barcode_seq/sequencing_data_Jan2023/'
    os.chdir(work_dir)

    # Set metatable directory
    meta_dir = 'metatable_Nov2022_all.txt'
    # meta_dir = 'metatable_Jan2023.txt'

    # Set output directory
    output_dir = 'Jan2023/'
    create_directories([output_dir])

    # Read the metatable that guides the column names to refer to
    df_meta = pd.read_csv(meta_dir, index_col='replicate', sep='\t')
    # All the column names are unique by default
    orf19_col = df_meta['orf19_column'].unique()[0]
    feature_col = df_meta['feature_column'].unique()[0]
    common_col = df_meta['common_column'].unique()[0]
    plate_col = df_meta['plate_column'].unique()[0]
    keep_col_list = [feature_col, common_col, plate_col]  # Column names to keep for some QC files
    plate_subset_col = df_meta['plate_subset'].unique()[0]  # Can choose a subset of plates from the read table
    if np.isnan(plate_subset_col):
        plate_subset_col = None
    else:
        plate_subset_col = plate_subset_col.split(';')  # Example entry: Plate 17;Plate18;Plate19

    groups = df_meta['group'].unique().tolist()
    for group in groups:
        print("Group:", group)
        group_dir = output_dir + group + '/'
        lfc_qc_dir = group_dir + 'LFC_QC/'
        dom_dir = group_dir + 'dominating_strains/'
        avgCorr_sumTop3_dir = group_dir + 'avgCorr_sumTop3/'
        mod_t_input_dir = group_dir + 'mod_t_input/'
        mod_t_output_dir = group_dir + 'mod_t_output/'
        create_directories([group_dir, lfc_qc_dir, dom_dir, avgCorr_sumTop3_dir, mod_t_input_dir, mod_t_output_dir])

        df_group = df_meta[df_meta['group'] == group]
        input_reads_dir = df_group['dataset'].unique()[0]
        print("Reading the input raw reads file:", input_reads_dir)
        # Read the raw reads for that group
        df_reads = pd.read_csv(input_reads_dir, index_col=orf19_col, sep='\t')

        inoc_up = df_group.loc['inoculum', 'UP']
        inoc_dn = df_group.loc['inoculum', 'DN']
        rep_up_list = df_group['UP'].drop('inoculum').tolist()
        rep_dn_list = df_group['DN'].drop('inoculum').tolist()
        rep_up_list_LFC = [x + '_LFC' for x in rep_up_list]
        rep_dn_list_LFC = [x + '_LFC' for x in rep_dn_list]
        rep_list_LFC = rep_up_list_LFC + rep_dn_list_LFC

        # 1. Check inoculum fraction
        print("1. Checking inoculum fraction...")
        plot_inoc_qc(df=df_reads, inoc_list=[inoc_up, inoc_dn], up_list=rep_up_list, dn_list=rep_dn_list,
                     plates=plate_subset_col, plate_col=plate_col, read_thre_good=read_cutoff,
                     read_thre_bad=dropout_cutoff,
                     inoc_prefix=group, output_direc=lfc_qc_dir, plt_shape=qc_plot_shape)

        # 2. Compile LFC files
        print("2. Compiling LFC files...\n")
        df_group_LFC = get_LFC_with_mask_singleInoc(df_reads, keep_col_list, inoc_up, inoc_dn, rep_up_list, rep_dn_list,
                                                    plate_col=plate_col, plate=plate_subset_col,
                                                    output_lfc_direc=lfc_qc_dir + group, mask_val=read_cutoff,
                                                    make_QC=True)

        # 3. Identify dominating strains and generate the corresponding summary table
        print("3. Generating summary table for dominating strains...\n")
        dom_up_dir = dom_dir + group + '_UP_Top3_dominating_strains.csv'
        dom_dn_dir = dom_dir + group + '_DN_Top3_dominating_strains.csv'
        dom_concat_dir = dom_dir + group + '_Top3_dominating_strains.csv'

        df_group_topN_UP = get_topN_dom_strains(df_reads.reset_index().set_index(feature_col),
                                                plate_list=df_reads[plate_col].unique(), mice_list=rep_up_list,
                                                output_direc=dom_up_dir, thre=0,
                                                plate_col=plate_col, orf19_col=orf19_col, common_col=common_col)

        df_group_topN_DN = get_topN_dom_strains(df_reads.reset_index().set_index(feature_col),
                                                plate_list=df_reads[plate_col].unique(), mice_list=rep_dn_list,
                                                output_direc=dom_dn_dir, thre=0,
                                                plate_col=plate_col, orf19_col=orf19_col, common_col=common_col)

        df_group_topN = pd.concat([df_group_topN_UP, df_group_topN_DN])
        df_group_topN.to_csv(dom_concat_dir, index=None)

        # 4. Calculate the average LFC correlation among replicates and plot its relationship with the sum of Top 3
        # dominating strains
        print(
            "4. Calculating the average PCC among the LFC values of replicates and plotting their relationship with the sum of Top 3...\n\n")
        avgCorr_sumTop3_output_dir = avgCorr_sumTop3_dir + group + '_avgCorr_vs_LFC_mask' + str(read_cutoff)
        get_avgCorr_vs_sumTop3(lfc_direc=lfc_qc_dir + group + '_LFC.csv',
                               dom_up_direc=dom_up_dir, dom_dn_direc=dom_dn_dir,
                               m_lfc_list_up=rep_up_list_LFC, m_lfc_list_dn=rep_dn_list_LFC,
                               output_direc=avgCorr_sumTop3_output_dir,
                               plate_col=plate_col, sort_cols=[plate_col, 'Mouse'])

    print("Step 1 of the pipeline finished.\nTo proceed, please follow the following steps:")
    print("  (1) Browse the QC plots and define the thresholds for average LFC correlation & dominating strains")
    print("  (2) Pick the thresholds, fill them into the metatable, and run the code for Step 2")
