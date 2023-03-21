# %%
import os
import pandas as pd 
import json
from MetaRes import get_Meta_res
from MetaRes2Rep import Res2Rep
from MetaRep2Excel import Rep2Dict, Rep2Excel


# %%
def get_micro(sample_name, p_file, g_file, s_file):
    m_dict = {}
    def get_level(file, level):
        df = pd.read_csv(file, sep='\t', index_col=0)
        df['tax_level'] = level
        df = df[[sample_name,'tax_level']]
        return df
    df_m = pd.DataFrame()
    for f, level in zip([p_file, g_file, s_file],['p','g','s']):
        df = get_level(f, level)
        df_m = pd.concat([df_m, df])

    m_dict['relative_abundance'] = df_m[sample_name].to_dict()
    m_dict['tax_level'] = df_m['tax_level'].to_dict()

    return m_dict

def get_alpha(sample_name, alpha_file):
    df = pd.read_csv(alpha_file, sep='\t', index_col=0)[sample_name]
    return df.to_dict()


def get_res_old(sample_name, p_file, g_file, s_file, alpha_file):
    total_res_dict = {}
    tax_value_dict = get_micro(sample_name, p_file, g_file, s_file)
    alpha_dict = get_alpha(sample_name, alpha_file)
    total_res_dict = {'tax_value_dict':tax_value_dict, 'alpha_dict':alpha_dict}
    return total_res_dict

# %%
def check_old_samples(config_file, db_dir, sample_id, latest_sample_n=3, outputdir='.'):
    p_file = os.path.join(db_dir, 'phylum')
    g_file = os.path.join(db_dir, 'genus')
    s_file = os.path.join(db_dir, 'sp')
    alpha_file = os.path.join(db_dir, 'alph')

    df_samples = pd.read_excel(os.path.join(db_dir, 'ID.xlsx'), index_col=0)
    sample_names = df_samples[df_samples['会员ID'] == sample_id].sort_values(by='采样时间', ascending=False).head(latest_sample_n)['样本ID'].to_list()

    rep_dict_list = []
    if sample_names:
        for sample_name in sample_names:
            total_res_dict = get_res_old(sample_name, p_file, g_file, s_file, alpha_file)
            tmp_sample_res_file = os.path.join(outputdir, f'{sample_name}.res.json')
            tmp_sample_rep_file = os.path.join(outputdir, f'{sample_name}.rep.json')
            with open(tmp_sample_res_file,'w') as h:
                json.dump(total_res_dict, h, ensure_ascii=False, indent=2)
            Res2Rep(tmp_sample_res_file, config_file, tmp_sample_rep_file)
            tmp_rep_dict = Rep2Dict(tmp_sample_res_file, tmp_sample_rep_file, config_file)
            rep_dict_list.append(tmp_rep_dict)
    return rep_dict_list

# %%
def check_samples(config_file, db_dir, sample_id, latest_sample_n=3, outputdir='.', old_db_dir=None):
    df_samples = pd.read_csv(os.path.join(db_dir, 'ID.csv'),index_col=0, header=None).drop_duplicates(subset=1)
    df_samples = df_samples.sort_values(1, ascending=False)
    sample_names = []
    if sample_id in df_samples.index:
        df_samples = df_samples.loc[sample_id]
        if df_samples.shape[0] > 1:
            sample_names = df_samples[1].head(latest_sample_n).to_list()
        else:
            sample_names = df_samples.to_list()
    rep_dict_list = []

    if sample_names:
        for sample_name in sample_names:
            sample_dir = os.path.join(db_dir, sample_name[:5], sample_name[:8], sample_name)
            rgi_res_file = os.path.join(sample_dir, 'card.txt')
            vf_res_file = os.path.join(sample_dir, 'vfdb.diamond.txt')
            raw2orf_file = os.path.join(sample_dir, 'raw2orf.bam.idxstats')
            metaphlan_res_file = os.path.join(sample_dir, 'final_metaphlan_bugs_list.tsv')
            pathabundance_res_file = os.path.join(sample_dir,'final_pathabundance_unstratified.tsv')
            tmp_sample_res_file = os.path.join(outputdir, f'{sample_name}.res.json')
            tmp_sample_rep_file = os.path.join(outputdir, f'{sample_name}.rep.json')
            get_Meta_res(tmp_sample_res_file, metaphlan_res_file, raw2orf_file, rgi_res_file, vf_res_file, pathabundance_res_file)
            Res2Rep(tmp_sample_res_file, config_file, tmp_sample_rep_file)
            tmp_rep_dict = Rep2Dict(tmp_sample_res_file, tmp_sample_rep_file, config_file)
            rep_dict_list.append(tmp_rep_dict)
    
    if len(sample_names) < latest_sample_n:
        old_res = latest_sample_n - len(sample_names)
        if old_db_dir:
            tmp_old_rep_dict = check_old_samples(config_file, old_db_dir, sample_id, latest_sample_n=old_res, outputdir=outputdir)
            rep_dict_list = rep_dict_list + tmp_old_rep_dict
    rep_dict_list.reverse()

    return rep_dict_list


# %%

# sample_id = 'nlwushubin'
# config_file = 'network_test_data.xlsx'
# db_dir = 'test/db/'
# old_db_dir = 'test/old_db/'
# rep_dict_list = check_samples(config_file, db_dir, sample_id, old_db_dir=old_db_dir)

# final_out_excel = 'test_out.xlsx'
# Rep2Excel(final_out_excel, rep_dict_list)

