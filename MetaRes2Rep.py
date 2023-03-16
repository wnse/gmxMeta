# %%
import os
import json
import pandas as pd
from FunctionNet.FunctionNet import FunctionNet

def get_sourc(excel_file):
    df_raw = pd.DataFrame()
    df_out = pd.DataFrame()
    path = os.path.realpath(os.path.split(excel_file)[0])
    reader = pd.ExcelFile(excel_file)
    if 'source' in reader.sheet_names:
        df_raw = pd.read_excel(reader, 'source')
    else:
        print(f"'source' not a sheet in {excel_file}")
        return None
    
    types = []
    if not df_raw.empty:
        if 'type' not in df_raw.columns:
            print(f"'type' not in source columns")
            return None
        if 'source' not in df_raw.columns:
            print(f"'source' not in source columns")
        else:
            types = df_raw['type'].unique()
        
    for type in types:
        df_tmp = df_raw.loc[df_raw['type']==type,:].copy()
        if type.upper() == 'VALUE':
            df_tmp['from'] = 'source'
            df_out = pd.concat([df_out,df_tmp])
        if type.upper() == 'SHEET':
            for i in df_tmp.index:
                source_tmp = df_tmp.loc[i,'source']
                if source_tmp in reader.sheet_names:
                    df_tmp_tmp = pd.read_excel(reader, source_tmp).astype(str)
                    df_tmp_tmp['from'] = source_tmp
                    if 'source' not in df_tmp_tmp.columns:
                        print(f"'source' not in {source_tmp} columns")
                    else:
                        df_out = pd.concat([df_out, df_tmp_tmp])
                else:
                    print(f"{source_tmp} not sheet in {excel_file}")
        if type.upper() == 'CSV':
            for i in df_tmp.index:
                source_tmp = df_tmp.loc[i, 'source']
                file_tmp = os.path.join(path, source_tmp)
                if os.path.isfile(file_tmp):
                    df_tmp_tmp = pd.read_csv(file_tmp, sep='\t').astype(str)
                    df_tmp_tmp['from'] = source_tmp
                    if 'source' not in df_tmp_tmp.columns:
                        print(f"'source' not in {file_tmp} columns")
                    else:
                        df_out = pd.concat([df_out, df_tmp_tmp])
                else:
                    print(f"{file_tmp} not csv in {path}")
        
    return df_out


def Res2Rep(res_json_file, config_file, out_rep_json=None):
    if not out_rep_json:
        out_rep_json=os.path.splitext(res_json_file)[0]+'.rep.json'
    reader = pd.ExcelFile(config_file)
    df_edge = pd.read_excel(reader,'edge').fillna('').astype(str)
    df_target = pd.read_excel(reader,'target').fillna('').astype(str)
    # df_source = pd.read_excel(reader, 'source')
    df_source = get_sourc(config_file)
    fn = FunctionNet(df_edge, df_target, df_source)

    with open(res_json_file,'rt') as h:
        total_res_dict_tmp = json.load(h)

    total_res_dict = {}
    for k, v in total_res_dict_tmp.items():
        if k == 'tax_value_dict':
            total_res_dict.update(v['relative_abundance'])
        else:
            total_res_dict.update(v)

    # fn.check()
    total_rep_dict = fn.compute_all_target(total_res_dict)
    with open(out_rep_json,'w') as h:
        json.dump(total_rep_dict, h, ensure_ascii=False, indent=2)

# %%
# file = 'network_test_data.xlsx'
# res_file = 'total_res_dict.json'
# Res2Rep(res_file, file)


