import networkx as nx
import pandas as pd
import json
import matplotlib.pyplot as plt
import pandas as pd
from FunctionNet.Function import Function

class FunctionNet(Function):
    # plt.rcParams['font.sans-serif'] = ['SimHei']
    def __init__ (self, df_edge, df_target=pd.DataFrame(), df_source=pd.DataFrame()):
        if self.check_repeat(df_edge, ['target','source'], 'edge'):
            self.g = nx.from_pandas_edgelist(df_edge, edge_attr=True)
            self.info = df_edge
        else:
            self.g = None
        if not df_target.empty:
            if 'target' in df_target.columns:
                if self.check_repeat(df_target, 'target', 'target'):
                    self.target = df_target.set_index('target').to_dict(orient='index')
                else:
                    self.target = {}
            else:
                print(f"'target' not in columns:{df_target.columns}")
        if not df_source.empty:
            if 'source' in df_source.columns:
                if self.check_repeat(df_source, 'source', 'source'):
                    self.source = df_source.set_index('source').to_dict(orient='index')
                else:
                    self.source = {}
            else:
                print(f"'source' not in columns:{df_source.columns}")

    def check_repeat(self, df, col, tag):
        repeat = df[df[col].duplicated()][col].values
        if len(repeat) > 0:
            for e in repeat:
                print(f"{tag} repeated: {e}")
            return None
        else:
            return 1

    def check(self):
        check_res = 'OK'
        for edge in self.g.edges():
            function_name = self.g.edges[edge].get('function')
            if function_name:
                try:
                    call_func = getattr(self, function_name)
                except Exception as e:
                    print(f"edge {edge} : {e}")
                    check_res = 'Not OK'
            else:
                print(f"'function' not in edge {edge} attributs")
                check_res = 'Not OK'
        for function_name in set([v['function'] for i,v in self.target.items()]):
            if function_name:
                try:
                    call_func = getattr(self, function_name)
                except Exception as e:
                    print(f"target function {function_name} : {e}")
                    check_res = 'Not OK'
            else:
                print(f"target function 【{function_name}】 not exists ")
                check_res = 'Not OK'
        for node in self.g.nodes():
            if node not in self.target.keys() and node not in self.source.keys():
                print(f"{node} not in tarrget and source")
                check_res = 'Not OK'
        
        print(check_res)

    def draw_networkx(self, outpng:str='FuntionNet.png'):
        plt.subplots(figsize=(10,10))
        nx.draw_networkx(self.g)
        plt.savefig(outpng, transparent=True, dpi=300)

    def list_all_source(self, target, source_data=None):
        def check_target_in_source(target):
            if target in self.source.keys():
                source_value = 0
                if source_data:
                    if target in source_data.keys():
                        print(f"{target}:{source_data[target]}")
                        source_value = source_data[target]
                    else:
                        print(f"without value in source data for {target}")
                return source_value
            else:
                return None
        check_target = check_target_in_source(target)
        if check_target == None:
            if target in self.target.keys():
                check_target_source_pass = []
                check_target_source_fail = []
                for layer in nx.bfs_successors(self.g, source=target):
                    tmp_target, tmp_source = layer
                    if tmp_target == target or (tmp_target in check_target_source_fail):
                        for ts in tmp_source:
                            check_source = check_target_in_source(ts)
                            if check_source == None:
                                check_target_source_fail.append(ts)
                            else:
                                check_target_source_pass.append(ts)
                                print(f"{ts}:{check_source}")
            else:
                print(f"{target} not in targets")
                return None
            return check_target_source_pass
        else:
            return check_target

    def compute(self, source_value, function_name, function_arg):
        try:
            call_func = getattr(self, function_name)
            if function_arg:
                args_list = function_arg.split('|')
                function_out = call_func(source_value,*args_list)
            else:
                function_out = call_func(source_value)
            return function_out
        except Exception as e:
            print(f"{function_name} {function_arg}: {e}")
            return None

    def compute_target(self, target, source_data):
        for_compute = self.info[self.info['target'] == target].set_index('source').to_dict(orient='index')
        compute_source_dict = {}
        # compute_target_dict = {}
        for s, info in for_compute.items():
            compute_source_dict[s] = info
            function_name = info.get('function')
            function_arg = info.get('arguments')
            if s in self.source.keys():
                # compute_source_dict[s]['source_result'] = source_data.get(s, 0)
                if function_name:
                    # if s not in source_data.keys():
                        # print(f"{s} not in source data, will be assigned 0")
                    compute_source_dict[s]['source_result'] = self.compute(source_data.get(s,0), function_name, function_arg)
                else:
                    print(f"function name error [{target} {s}]: {function_name}")
            elif s in self.target.keys():
            # else:
                # compute_source_dict[s]['target_result'] = compute_source_dict[s].get('target_result', None)
                _, compute_source_dict[s]['target_result'] = self.compute_target(s, source_data)
                compute_source_dict[s]['source_result'] = self.compute(compute_source_dict[s]['target_result'], function_name, function_arg)


        if target in self.target.keys():
            function_name = self.target[target].get('function')
            function_arg = self.target[target].get('arguments')
            if function_name:
                # print(f'comput target {target} {function_name}')
                target_result =  self.compute(compute_source_dict, function_name, function_arg)
            else:
                print(f"function name error [{target}]: {function_name}")

        return compute_source_dict, target_result

    def compute_all_target(self, source_data):
        # pass
        computed_result = {}
        for t in self.target.keys():
            if t in computed_result.keys():
                next
            else:
                compute_source, target_result  = self.compute_target(t, source_data)
                computed_result[t] = target_result
                for s_key, s_info in compute_source.items():
                    if s_key in self.target.keys() and (s_key not in computed_result.keys()):
                        computed_result[s_key] = s_info['target_result']
                        # print(f"{s_key}:{s_info['target_result']}")
        return computed_result


