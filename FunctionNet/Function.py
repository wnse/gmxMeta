from typing import List, Dict
import numpy as np

class Function(object):
    def __init__(self):
        pass

    # def RawRisk(self, x:float, range_min:float, range_max:float) -> int:
    #     res = 0
    #     if x <= float(range_min):
    #         res = -1
    #     if x > float(range_max):
    #         res = 1
    #     return res

    def RawRisk(self, x:float, *args) -> int:
        range_list = np.array(args, dtype=float)
        res_list = [-1] + [0]*(len(range_list)-1) +[1]
        return res_list[np.searchsorted(sorted(range_list), x)]
        
    def RawValue(self, x:float) -> float:
        return x
    
    def Sum(self, source:Dict) -> float:
        function_result = 0
        for s, info in source.items():
            function_result += info['source_result']
        return function_result

    # def OneForAll(self, source:Dict, cutoff, risk=1) -> int:
    #     function_result = 0
    #     for s, info in source.items():
    #         if info['source_result'] != float(cutoff):
    #             function_result = risk
    #     return function_result
    
    def OneForAll(self, source:Dict, cutoff, risk=1) -> int:
        function_result = {}
        function_result['source_data'] = {}
        tmp_res = 0
        for s, info in source.items():
            function_result['source_data'][s] = info['source_result']
            if info['source_result'] != float(cutoff):
                tmp_res = risk
        function_result['target_result'] = tmp_res
        return function_result
    
    def SumForAll(self, source:Dict, range_min:float, range_max:float) -> int:
        function_result = 0
        for s, info in source.items():
            function_result += info['source_result']
        function_result = self.RawRisk(function_result, range_min, range_max)
        return function_result
    
    def MaxforAll(self, source:Dict) -> int:
        function_result_val = 0
        function_result = ''
        for s, info in source.items():
            if info['source_result'] >= function_result_val:
                function_result_val = info['source_result']
                function_result = s
        return function_result
    
    def OneForOne(self, source:Dict) -> dict:
        tmpi = [i for i, v in source.items()]
        tmpv = [v['source_result'] for i, v in source.items()]
        function_result = dict(zip(tmpi, tmpv))
        return function_result
    


