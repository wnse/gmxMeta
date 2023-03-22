# %%
import pandas as pd
import json
import os
from MetaAlphaPng import alph_png

def read_res(res_json_file, rep_json_file, config_excel):
    rep_json = {}
    with open(rep_json_file, 'rt') as h:
        rep_json = json.load(h)

    res_json = {}
    with open(res_json_file, 'rt') as h:
        res_json = json.load(h)

    df_micro_meta = pd.read_excel(config_excel, 'Microbe')
    df_edge = pd.read_excel(config_excel, 'edge')
    df_target = pd.read_excel(config_excel, 'target')

    return res_json, rep_json, df_micro_meta, df_edge, df_target

def get_micro_res(tag, res_json, rep_json, df_micro_meta):
    df_rep = pd.DataFrame.from_dict(rep_json[tag], orient='index', columns=['res'])
    df_rep = pd.merge(df_rep, df_micro_meta,left_index=True, right_on='source', how='left')
    df_rep['risk'] = df_rep['res'].map({0:'正常',1:'高',-1:'低'})
    df_rep['result'] = df_rep['source'].map(res_json['tax_value_dict']['relative_abundance'])
    df_rep['result'] = df_rep['result'].fillna(0).round(3).astype(str) + '%'
    df_rep = df_rep.rename(columns={'拉丁文名':tag}).fillna('')

    micro_res = '、'.join(df_rep[df_rep['res']!=0][tag].to_list())
    df_rep = df_rep.drop(['res','source'],axis=1)
    return micro_res, df_rep

def get_micro_res_all(res_json, rep_json, df_micro_meta, micro_tags):
    micro_df = []
    micro_res = []
    for tag in micro_tags:
        res_tmp, df_tmp =  get_micro_res(tag, res_json, rep_json, df_micro_meta)
        micro_df.append(df_tmp.set_index(tag))
        micro_res.append(res_tmp)
    return micro_res, micro_df


def get_enterotype_res(res_json, rep_json, df_micro_meta, df_edge):
    enter_source = df_edge[df_edge['target'] == '肠型'].set_index('source')
    enterotype_res = enter_source.loc[rep_json['肠型']]['tag']
    df_enterotype_res = df_micro_meta[df_micro_meta['source'].isin(enter_source.index)].set_index('source').rename_axis('').copy()
    df_enterotype_res['检测结果（%）'] = df_enterotype_res.index.map(res_json['tax_value_dict']['relative_abundance'])
    if '中文名称' in df_enterotype_res.columns:
        df_enterotype_res = df_enterotype_res[['中文名称','检测结果（%）']]
    df_enterotype_res = df_enterotype_res.fillna(0)
    return enterotype_res, df_enterotype_res

def get_diversity_res(res_json, rep_json, df_edge):
    df_diversity = df_edge[df_edge['target'] == 'α多样性'].copy()

    diversity_res = pd.DataFrame.from_dict(rep_json['α多样性'],orient='index',columns=['risk'])
    diversity_res['tag'] = diversity_res.index.map(df_diversity.set_index('source')['tag'].to_dict())
    risk_exp = {0:'处于中国健康人群范围内',-1:'低于中国健康人群范围',1:'处于中国健康人群范围内'}
    diversity_res['risk_exp'] = diversity_res['risk'].map(risk_exp)
    if len(diversity_res['risk_exp'].unique())==1:
        diversity_res_exp = '、'.join(diversity_res['tag'].to_list()) + diversity_res['risk_exp'].unique()[0]
    else:
        diversity_res_tmp = diversity_res[diversity_res['risk'] == -1]
        diversity_res_exp = '、'.join(diversity_res_tmp[['tag','risk_exp']].apply(''.join,axis=1))

    df_diversity['参考值'] = df_diversity['arguments'].str.split('|').str[0]
    df_diversity['result'] = df_diversity['source'].map(res_json['alpha_dict'])
    df_diversity['result'] = df_diversity['result'].round(2)
    df_diversity = df_diversity.rename(columns={'tag':'α多样性'})
    df_diversity = df_diversity[['α多样性','参考值','result']].set_index('α多样性')

    return diversity_res_exp, df_diversity

def get_disease_res(res_json, rep_json, df_micro_meta, df_edge, df_target):
    df_disease = df_target[df_target['tag'].isin(['菌群与疾病风险评估'])].copy()
    for i in df_disease['target'].unique():
        df_disease.loc[df_disease['target']==i, 'risk'] = rep_json[i]['target_result']
    # df_disease['risk'] = df_disease['target'].map(rep_json)
    risk_exp = {1:'高度风险',0:'低风险'}
    df_disease['risk_exp'] = df_disease['risk'].map(risk_exp)
    tmp_exp = ['高度风险', '低风险']
    disease_res_exp = []
    for tmp in tmp_exp:
        disease_res_exp.append("、".join(df_disease.loc[df_disease['risk_exp'] == tmp,'target'].to_list()))

    df_disease_res = df_edge[df_edge['target'].isin(df_disease['target'])].copy()
    df_disease_res['result'] = df_disease_res['source'].map(res_json['tax_value_dict']['relative_abundance']).fillna(0)
    df_disease_res['检测含量'] = df_disease_res['result']
    for i in df_disease_res['target'].unique():
        df_disease_res.loc[df_disease_res['target']==i, 'risk'] = df_disease_res.loc[df_disease_res['target']==i, 'source'].map(rep_json[i]['source_data'])
    
    mask = df_disease_res['source'].isin(['短链脂肪酸合成'])
    df_disease_res.loc[~mask, '检测含量'] = df_disease_res.loc[~mask, '检测含量'].round(2).astype(str) + '%'
    risk_exp = {1:'偏高',0:'正常',-1:'偏低'}
    df_disease_res['结果'] = df_disease_res['risk'].map(risk_exp)

    df_disease_res['微生物种类'] = df_disease_res['source'].map(df_micro_meta.set_index('source')['拉丁文名'])
    check_na = df_disease_res['微生物种类'].isna()
    df_disease_res.loc[check_na,'微生物种类'] = df_disease_res.loc[check_na, 'source']
    df_disease_res = df_disease_res.rename(columns={'tag':'健康人群含量'})
    df_disease_res = df_disease_res[['微生物种类','检测含量','健康人群含量','结果']].set_index('微生物种类')

    return disease_res_exp, df_disease_res

def get_vf_res(rep_json, df_target):
    df_vf = df_target[df_target['tag'].isin(['肠道微生物组毒力因子分析'])].copy()
    df_vf['risk'] = df_vf['target'].map(rep_json)
    risk_exp = {1:'检出',0:'未检出'}
    df_vf['检测结果'] = df_vf['risk'].map(risk_exp)
    df_vf['毒力因子'] = df_vf['target']

    vf_list = '、'.join(df_vf[df_vf['risk'] == 1]['target'].to_list())
    return vf_list, df_vf[['毒力因子','检测结果']].set_index('毒力因子')

def get_rg_res(rep_json, df_target):
    df_rg = df_target[df_target['tag'].isin(['肠道微生物组抗生素抗性基因分析'])].copy()
    df_rg['risk'] = df_rg['target'].map(rep_json)
    risk_exp = {1:'检出',0:'未检出'}
    df_rg['检测结果'] = df_rg['risk'].map(risk_exp)
    df_rg['抗生素类别'] = df_rg['target']

    rg_list = '、'.join(df_rg[df_rg['risk'] == 1]['target'].to_list())

    return rg_list, df_rg[['抗生素类别','检测结果']].set_index('抗生素类别')


def get_pathway_res(rep_json, df_target):
    df_pathway = df_target[df_target['tag'].isin(['代谢能力评估'])].copy()
    df_pathway['risk'] = df_pathway['target'].map(rep_json)
    risk_exp = {1:'偏高',0:'正常',-1:"偏低"}
    df_pathway['检测结果'] = df_pathway['risk'].map(risk_exp)
    dict_pathway = df_pathway.set_index('target')['检测结果'].to_dict()

    return dict_pathway



def get_tax_top(res_json, top=10):
    df_tax = pd.DataFrame.from_dict(res_json['tax_value_dict']['tax_level'], orient='index', columns=['tax_level'])
    df_ratio = pd.DataFrame.from_dict(res_json['tax_value_dict']['relative_abundance'], orient='index', columns=['relative_abundance'])
    df_tax = pd.merge(df_tax, df_ratio, left_index=True, right_index=True)
    
    tax_levels = ['p','g','s']
    top_tax_list = []
    for tax_level in tax_levels:
        df_tax_tmp = df_tax[df_tax['tax_level'] == tax_level]['relative_abundance'].sort_values(ascending=False)
        df_top = df_tax_tmp.head(top)
        rest = sum(df_tax_tmp.to_list()[10:])
        if rest:
            df_top.loc['Rest'] = rest
        df_top = df_top.rename_axis('clade_name')
        top_tax_list.append(df_top)
    
    return top_tax_list


# %%

def Rep2Dict(res_json_file, rep_json_file, config_file):
    config_excel = pd.ExcelFile(config_file)
    res_json, rep_json, df_micro_meta, df_edge, df_target = read_res(res_json_file, rep_json_file, config_excel)
    ### 各类微生物
    micro_tags = ['核心菌属','有害菌属', '其他菌属', '有益微生物', '有害微生物', '条件致病菌', '其他微生物']
    micro_res_list, micro_res_df_list = get_micro_res_all(res_json, rep_json, df_micro_meta, micro_tags)
    top_tax_res_df_list = get_tax_top(res_json)
    ### 肠型
    enterotype_res, enterotype_res_df = get_enterotype_res(res_json, rep_json, df_micro_meta, df_edge)
    ### 多样性
    diversity_res, diversity_res_df = get_diversity_res(res_json, rep_json, df_edge)
    ### 疾病风险
    disease_res_list, disease_res_df = get_disease_res(res_json, rep_json, df_micro_meta, df_edge, df_target)
    ### 代谢通路
    pathway_res_dict = get_pathway_res(rep_json, df_target)
    ### 毒力因子
    vf_res, vf_res_df = get_vf_res(rep_json, df_target)
    ### 抗性基因
    rg_res, rg_res_df = get_rg_res(rep_json, df_target)
    ### summary
    df_summary = {}
    df_summary['菌群结构'] = '门水平，属水平，种水平'
    for tag, res in zip(micro_tags, micro_res_list):
        df_summary[tag] = res
    df_summary['肠型'] = enterotype_res
    df_summary['微生物组均匀度'] = diversity_res
    df_summary['高度风险'] = disease_res_list[0]
    df_summary['低风险'] = disease_res_list[1]
    pathway_tags = ['能量代谢', '碳水化合物代谢', '脂类代谢', '蛋白质代谢', '必须氨基酸合成', '维生素合成', 
                    '短链脂肪酸合成', '抗氧化能力', '有害物代谢', '脂多糖合成']
    for tag, res in zip(pathway_tags, pathway_res_dict):
        if tag in ['能量代谢', '碳水化合物代谢', '脂类代谢', '蛋白质代谢']:
            df_summary[tag] = '正常'
        else:
            df_summary[tag] = pathway_res_dict[tag]
    df_summary['抗生素暴露史'] = rg_res
    df_summary['毒力因子'] = vf_res
    df_summary['总结'] = ''
    df_summary = pd.DataFrame.from_dict(df_summary, orient='index')
    ### 输出dict
    out_json_dict = {}
    out_df_list = [df_summary, enterotype_res_df, diversity_res_df, rg_res_df, vf_res_df, disease_res_df,
                    top_tax_res_df_list[0], top_tax_res_df_list[1], top_tax_res_df_list[2]]
    out_tag_list = ['summary', '肠型', '多样性', '抗生素抗性基因', '毒力因子', '疾病风险', 'Phylum','Genus','Species']
    for df, tag in zip(out_df_list, out_tag_list):
        out_json_dict[tag] = df.reset_index().drop_duplicates(subset=df.reset_index().columns).to_dict()
    for i, tag in enumerate(micro_tags):
        out_json_dict[tag] = micro_res_df_list[i].reset_index().to_dict()

    return out_json_dict


def Rep2Excel(out_excel, res_dict_list, min_sample=3):
    if not res_dict_list:
        return None
    # out_excel = 'total_rep.xlsx'
    micro_tags = ['核心菌属','有害菌属', '其他菌属', '有益微生物', '有害微生物', '条件致病菌', '其他微生物']
    ### 获取最后一次结果
    res_dict = res_dict_list[-1]
    in_df_list = []
    in_tag_list = ['summary', '肠型', '多样性', '抗生素抗性基因', '毒力因子', '疾病风险']
    for i, tag in enumerate(in_tag_list):
        in_df_list.append(pd.DataFrame.from_dict(res_dict[tag]))
    df_summary, enterotype_res_df, diversity_res_df, rg_res_df, vf_res_df, disease_res_df = in_df_list
    top_tax_res_df_list = []
    for level in ['Phylum','Genus','Species']:
        top_tax_res_df_list.append(pd.DataFrame.from_dict(res_dict[level]))
    ### 整理多样性多次结果
    for i, res_dict_tmp in enumerate(res_dict_list):
        df_tmp = pd.DataFrame.from_dict(res_dict_tmp['多样性'])
        diversity_res_df.insert(i+1, f"T{i+1}", df_tmp['result'])
    diversity_res_df = diversity_res_df.drop('result', axis=1)
    ### 补齐min_sample次
    if len(res_dict_list) < min_sample:
        for i in range(len(res_dict_list), min_sample):
            diversity_res_df.insert(i+1, f"T{i+1}", '')

    ### 整理各种类微生物多次结果
    micro_res_df_list = []
    for i, tag in enumerate(micro_tags):
        df_tmp_micro = pd.DataFrame.from_dict(res_dict[tag])
        for j, res_dict_tmp in enumerate(res_dict_list):
            df_tmp = pd.DataFrame.from_dict(res_dict_tmp[tag])
            df_tmp_micro.insert((j+1)*2, f"T{j+1}", df_tmp['result'])
            df_tmp_micro.insert((j+1)*2+1, f"T{j+1}结果", df_tmp['risk'])
        if len(res_dict_list) < min_sample:
            for j in range(len(res_dict_list), min_sample):
                df_tmp_micro.insert((j+1)*2, f"T{j+1}", '')
                df_tmp_micro.insert((j+1)*2+1, f"T{j+1}结果", '')
        df_tmp_micro = df_tmp_micro.drop(['result','risk'], axis=1)
        micro_res_df_list.append(df_tmp_micro)
    ### 输出excel
    try:
        writer = pd.ExcelWriter(out_excel)
        df_summary.T.to_excel(writer, '检测结果汇总', index=False, header=False)
        for i, tag in enumerate(micro_tags):
            micro_res_df_list[i].to_excel(writer, tag, index=False)
        enterotype_res_df.to_excel(writer, '肠型', index=False)
        diversity_res_df.to_excel(writer, '多样性', index=False)
        rg_res_df.to_excel(writer, '抗生素抗性基因', index=False)
        vf_res_df.to_excel(writer, '毒力因子', index=False)
        disease_res_df.to_excel(writer, '疾病风险', index=False)
        for level, tax in zip(['Phylum','Genus','Species'], top_tax_res_df_list):
            tax.to_excel(writer, level, index=False)
        writer.close()
    except Exception as e:
        print(f"ERROR write excel {out_excel} {e}")

    try:
        # out_png = os.path.splitext(out_excel)[0]+'.png'
        out_png = os.path.join(os.path.split(out_excel)[0], '4.alph.png')
        alph_png(diversity_res_df, out_png)
    except Exception as e:
        print(f"ERROR write png {out_png} {e}")
        


# %%
# res_json_file = 'test_res.json'
# rep_json_file = 'test_res.rep.json'
# config_file = 'network_test_data.xlsx'
# out_json_dict = Rep2Dict(res_json_file, rep_json_file, config_file)

# Rep2Excel('test_out.xlsx', [out_json_dict])

# %%
# with open('total_rep.json', 'w') as h:
    # json.dump(out_json_dict, h, ensure_ascii=False, indent=2)


