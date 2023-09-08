# Author: Xiang Zhang
# Contact: zhan6668@umn.edu
# Updated on 09/06/2023

# Step 2: After updating the metatable with the thresholds for average LFC correlation and dominating strains,
# we proceed with this step to generate the required input format for applying moderated t-test in R
# The user needs to manually set the directories for working and metatable.


from utils import *


if __name__ == '__main__':
    # Load the configuration from the JSON file
    with open('config.json', 'r') as config_file:
        config = json.load(config_file)

    # Access the parameters from the config dictionary
    read_cutoff = config['read_cutoff']  # Default is 50

    # Set working directory
    work_dir = config['work_dir']
    os.chdir(work_dir)

    # Set metatable directory
    meta_dir = config['meta_dir']

    # Set output directory
    output_dir = config['output_dir']

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
        #create_directories([group_dir, lfc_qc_dir, dom_dir, avgCorr_sumTop3_dir, mod_t_input_dir, mod_t_output_dir])

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

        dom_concat_dir = dom_dir + group + '_Top3_dominating_strains.csv'
        avgCorr_sumTop3_output_dir = avgCorr_sumTop3_dir + group + '_avgCorr_vs_LFC_mask' + str(read_cutoff)
        # Directory for the variance-normalized LFC file
        lfc_rm_norm_output_dir = lfc_qc_dir + group + '_LFC_rmTop3_'+str(dom_cutoff)+'_pcc'+str(corr_cutoff)+'_var_norm.txt'
        # Directory for the variance-normalized LFC values formatted into the R function input
        lfc_rm_norm_mod_output_dir = mod_t_input_dir + group + '_LFC_var_modt_input'

        print("Generating the finalized LFC files...")
        print("Here are the replicates excluded based on our selection of thresholds:")
        get_var_norm_mod_t_input_with_filter(dataset_type=group+'_',
                                             lfc_direc=lfc_qc_dir+group+'_LFC.csv',
                                             dom_direc=dom_concat_dir,
                                             avg_up_direc=avgCorr_sumTop3_output_dir+'_UP.txt',
                                             avg_dn_direc=avgCorr_sumTop3_output_dir+'_DN.txt',
                                             m_lfc_list_up=rep_up_list_LFC, m_lfc_list_dn=rep_dn_list_LFC,
                                             lfc_rm_norm_output_direc=lfc_rm_norm_output_dir,
                                             lfc_rm_norm_mod_output_direc=lfc_rm_norm_mod_output_dir,
                                             plate_col=plate_col, dom_thre=float(dom_cutoff), corr_thre=float(corr_cutoff),
                                             up_tag=config['up_tag'], dn_tag=config['dn_tag'])
        print("\n\n")

    print("Step 2 of the pipeline finished.\nTo proceed, please follow the following steps:")
    print("  (3) Run the code for Step 3 in R specifying the input directories to score the data")
    print("  (4) After getting the scores, run the code for Step 4 for post-processing and harvest the outcome")

