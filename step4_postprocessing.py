# Author: Xiang Zhang
# Contact: zhan6668@umn.edu
# Updated on 12/28/2022

# Step 4: Generate the summary files for the output from moderated t-test in R
# The user needs to manually set the directories for working, metatable and description files.


from utils import *


if __name__ == '__main__':
    # Set working directory
    work_dir = '/Users/zhangxiang/Documents/Research_CB/GeneEssentiality_CandidaAlbicans/barcode_seq/sequencingforseptembernoblesamples/'
    os.chdir(work_dir)

    # Set metatable directory
    meta_dir = 'metatable_Nov2022_all.txt'

    # Set output directory
    output_dir = 'Nov2022/'

    # Set description reference files:
    maporf19_direc = '/Users/zhangxiang/Documents/Research_CB/GeneEssentiality_CandidaAlbicans/barcode_seq/map_gene/ORF19_Assembly22_mapping.tab.txt'
    descri_direc = '/Users/zhangxiang/Documents/Research_CB/GeneEssentiality_CandidaAlbicans/barcode_seq/map_gene/C_albicans_SC5314_A22_current_chromosomal_feature.tab.txt'
    Scer_ortho_direc = '/Users/zhangxiang/Documents/Research_CB/GeneEssentiality_CandidaAlbicans/barcode_seq/map_gene/C_albicans_SC5314_S_cerevisiae_orthologs.txt'

    # Read the metatable that guides the column names to refer to
    df_meta = pd.read_csv(meta_dir, index_col='replicate', sep='\t')
    # All the column names are unique by default
    plate_col = df_meta['plate_column'].unique()[0]

    groups = df_meta['group'].unique().tolist()
    for group in groups:
        print("Group:", group)
        group_dir = output_dir + group + '/'
        lfc_qc_dir = group_dir + 'LFC_QC/'
        dom_dir = group_dir + 'dominating_strains/'
        avgCorr_sumTop3_dir = group_dir + 'avgCorr_sumTop3/'
        mod_t_input_dir = group_dir + 'mod_t_input/'
        mod_t_output_dir = group_dir + 'mod_t_output/'
        # create_directories([group_dir, lfc_qc_dir, dom_dir, avgCorr_sumTop3_dir, mod_t_input_dir, mod_t_output_dir])

        df_group = df_meta[df_meta['group'] == group]
        inoc_up = df_group.loc['inoculum', 'UP']
        inoc_dn = df_group.loc['inoculum', 'DN']
        rep_up_list = df_group['UP'].drop('inoculum').tolist()
        rep_dn_list = df_group['DN'].drop('inoculum').tolist()
        rep_up_list_LFC = [x + '_LFC' for x in rep_up_list]
        rep_dn_list_LFC = [x + '_LFC' for x in rep_dn_list]
        rep_list_LFC = rep_up_list_LFC + rep_dn_list_LFC

        dom_cutoff = df_group['dominating_cutoff'].unique()[0]
        corr_cutoff = df_group['correlation_cutoff'].unique()[0]

        # Directory for the variance-normalized LFC file
        lfc_rm_norm_output_dir = lfc_qc_dir + group + '_LFC_rmTop3_' + str(dom_cutoff) + '_pcc' + str(
            corr_cutoff) + '_var_norm.txt'

        print("Formatting the output files from moderated t-test...")
        df_modt_output = generate_mod_t_integrated_table(
            lfc_direc=lfc_rm_norm_output_dir,
            m_lfc_list=rep_list_LFC,
            mod_output_direc=mod_t_output_dir+group+'_LFC_var_modt_output/',
            output_direc=mod_t_output_dir+group+'_results.txt', plate_col=plate_col,
            df_maporf19_direc=maporf19_direc, df_descri_direc=descri_direc, df_Scer_ortho_direc=Scer_ortho_direc)

        print("\n")

    result_dir = output_dir + 'Summary_mod_t_results.xlsx'
    print("Generating the final summary Excel file", result_dir, "...\n")
    with pd.ExcelWriter(result_dir) as writer:
        for group in groups:
            group_dir = output_dir + group + '/'
            mod_t_output_dir = group_dir + 'mod_t_output/'
            df_temp = pd.read_csv(mod_t_output_dir+group+'_results.txt', sep='\t')
            df_temp.to_excel(writer, sheet_name=group, index=None)

    print("Step 4 of the pipeline finished.")
    print("Please refer to the summary Excel file under your working directory for the final results.")

