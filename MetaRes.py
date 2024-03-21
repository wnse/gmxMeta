# %%
import pandas as pd
import numpy as np
import json
import logging
import os

# %%
def cal_simpson(df_ratio, N=1000):
    return 1- (df_ratio*N * (df_ratio*N - 1)).sum()/(N*100*(N*100-1))

def cal_shannon(df_ratio):
    return -(df_ratio/100 * np.log2(df_ratio/100)).sum()

def get_tax_ratio(metaphlan_res_file):
    df_metaphlan_res = pd.read_csv(metaphlan_res_file, sep='\t', skiprows=4)
    df_metaphlan_res['taxonomy_level'] = df_metaphlan_res['#clade_name'].str.split('|').str[-1]
    df_metaphlan_res['tax_level'] = df_metaphlan_res['#clade_name'].str.split('|').str[-1].str.split('__').str[0]
    df_metaphlan_res['taxonomy'] = df_metaphlan_res['#clade_name'].str.split('|').str[-1].str.split('__').str[1]
    df_ratio = df_metaphlan_res[df_metaphlan_res['tax_level'] == 's']['relative_abundance']
    alpha_simpson = cal_simpson(df_ratio, N=1e6)    
    alpha_shannon = cal_shannon(df_ratio)

    taxonomy_dict = df_metaphlan_res.set_index('taxonomy')[['relative_abundance','tax_level']].to_dict()
    taxonomy_dict.update({'relative_abundance_level':df_metaphlan_res.set_index('taxonomy_level')['relative_abundance'].to_dict()})
    alpha_dict = {'Shannon':alpha_shannon, 'Simpson':alpha_simpson,}

    return taxonomy_dict, alpha_dict

def get_pathabundance_ratio(pathabundance_res_file, sample_name=None):
    if not sample_name:
        sample_name = 'relative_abundance'
    df_pathway_res = pd.read_csv(pathabundance_res_file, sep='\t', index_col=0)
    df_pathway_res.columns = [sample_name] + df_pathway_res.columns.tolist()[1:]
    df_pathway_res = (df_pathway_res/df_pathway_res.sum())[sample_name]
    return df_pathway_res.to_dict()

def get_orf_ratio(raw2orf_file):
    df_raw2orf = pd.read_csv(raw2orf_file, sep='\t', header=None)
    df_raw2orf.columns = ['orf','length','mapReads','unmapReads']
    df_raw2orf = df_raw2orf.set_index('orf')
    rawReads = df_raw2orf.sum()[['mapReads','unmapReads']].sum()
    df_orfRatio = df_raw2orf.drop('*')['mapReads'] / rawReads
    return df_orfRatio

def get_aro_ratio(raw2orf_file, rgi_res_file, sample_name=None):
    df_orfRatio = get_orf_ratio(raw2orf_file)
    df_rgi_res = pd.read_csv(rgi_res_file, sep='\t')
    df_rgi_res['orf'] = df_rgi_res['ORF_ID'].str.split(' ').str[0]
    df_rgi_res = df_rgi_res.set_index('orf')['ARO'].astype(str)
    df_aroRatio = pd.merge(df_rgi_res, df_orfRatio, left_index=True, right_index=True, how='left')
    df_aroRatio = df_aroRatio.groupby('ARO')['mapReads'].sum()
    df_aroRatio = df_aroRatio.drop(df_aroRatio[df_aroRatio==0].index)
    if sample_name:
        df_aroRatio = df_aroRatio.rename(sample_name)
    return df_aroRatio.to_dict()

def get_vf_ratio(raw2orf_file, vf_res_file, sample_name=None):
    df_orfRatio = get_orf_ratio(raw2orf_file)
    df_vf_res = pd.read_csv(vf_res_file, sep='\t', header=None)
    df_vf_res = df_vf_res.sort_values(by=[10,2,11], ascending=[True, False, False])
    df_vf_res = df_vf_res.drop_duplicates(subset=[0],keep='first').set_index(0)[1].rename('VF')
    df_vfRatio = pd.merge(df_vf_res, df_orfRatio, left_index=True, right_index=True, how='left')
    df_vfRatio = df_vfRatio.groupby('VF')['mapReads'].sum()
    df_vfRatio = df_vfRatio.drop(df_vfRatio[df_vfRatio==0].index)
    if sample_name:
        df_vfRatio = df_vfRatio.rename(sample_name)
    return df_vfRatio.to_dict()





# %%
def get_Meta_res(out_file, metaphlan_res_file=None, raw2orf_file=None, rgi_res_file=None, vf_res_file=None, pathabundance_res_file=None, sample_name=None):
    tax_value_dict = {}
    alpha_dict = {}
    if metaphlan_res_file:
        if os.path.isfile(metaphlan_res_file):
            tax_value_dict, alpha_dict = get_tax_ratio(metaphlan_res_file)
        else:
            logging.error(f"not exists {metaphlan_res_file}")

    rgi_res_dict = {}
    if raw2orf_file and rgi_res_file:
        if os.path.isfile(raw2orf_file):
            if os.path.isfile(rgi_res_file):
                rgi_res_dict = get_aro_ratio(raw2orf_file, rgi_res_file, sample_name=sample_name)
            else:
                logging.error(f"not exists {rgi_res_file}")
        else:
            logging.error(f"not exists {raw2orf_file}")

    vf_res_dict = {}
    if raw2orf_file and vf_res_file:
        if os.path.isfile(raw2orf_file):
            if os.path.isfile(vf_res_file):
                vf_res_dict = get_vf_ratio(raw2orf_file, vf_res_file, sample_name=sample_name)
            else:
                logging.error(f"not exists {vf_res_file}")
        else:
            logging.error(f"not exists {raw2orf_file}")

    pathway_res_dict = {}
    if pathabundance_res_file:
        if os.path.isfile(pathabundance_res_file):
            pathway_res_dict = get_pathabundance_ratio(pathabundance_res_file, sample_name=sample_name)
        else:
            logging.error(f"not exists {pathabundance_res_file}")
            
    total_res_dict = {"tax_value_dict":tax_value_dict,
                    "alpha_dict":alpha_dict,
                    "rgi_res_dict":rgi_res_dict, 
                    "vf_res_dict":vf_res_dict, 
                    "pathway_res_dict":pathway_res_dict}
    with open(out_file,'w') as h:
        json.dump(total_res_dict, h, ensure_ascii=False, indent=2)

# %%

# sample_name = 'ma-mix-2159yk'
# rgi_res_file = '/mnt/c/work/metagenome/output/ma-mix-2159yk_out/card.txt'
# vf_res_file = '/mnt/c/work/metagenome/output/ma-mix-2159yk_out/vfdb.diamond.txt'
# raw2orf_file = '/mnt/c/work/metagenome/output/ma-mix-2159yk_out/raw2orf.bam.idxstats'

# metaphlan_res_file = ('/mnt/c/work/metagenome/output/ma-mix-2159yk_out/'
#                       'ma-mix-2159yk_FKDL190725375-1a-9_1_humann_temp/ma-mix-2159yk_FKDL190725375-1a-9_1_metaphlan_bugs_list.tsv')
# pathabundance_res_file = ('/mnt/c/work/metagenome/output/ma-mix-2159yk_out/'
#                           'ma-mix-2159yk_FKDL190725375-1a-9_1_pathabundance_unstratified.tsv')

# get_Meta_res('total_res.json', metaphlan_res_file, raw2orf_file, rgi_res_file, vf_res_file, pathabundance_res_file)


