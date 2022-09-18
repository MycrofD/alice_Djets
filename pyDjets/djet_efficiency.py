import contextlib
from pprint import pprint
import json
from os.path import exists

import uproot
import numpy as np
from central_config import f_dmeson_species, NDMC
from run_settings import run_mc_eff_file, run_is_prefix, run_is_postfix


def djet_efficiency(
    mc_eff_file: str = run_mc_eff_file,
    is_prefix: bool = run_is_prefix,
    is_postfix: bool = run_is_postfix,
    fix_name: str = "",
):
    root_file = uproot.open(mc_eff_file)
    directory_file = root_file["DmesonsForJetCorrelations"]

    hist_list, sparse_dphiz = [], []
    for i in range(NDMC):
        if not is_prefix and not is_postfix:
            hist_list.append(directory_file[f"histos{f_dmeson_species}MBN{i}MCrec"])
        if not is_prefix and is_postfix:
            hist_list.append(
                directory_file[f"histos{f_dmeson_species}MBN{i}{fix_name}MCrec"]
            )
        if is_prefix and is_postfix:
            hist_list.append(
                directory_file[f"histos{f_dmeson_species}{fix_name}MBN{i}MCrec"]
            )
        if is_prefix and not is_postfix:
            raise "Postfix has to be true if prefix is!"

        # print(vars(hist_list[i]))
        # print('dir,')
        # pprint(dir(hist_list[0]))
        # sparse_dphiz.append(np.array(hist_list[i].bases))
        # print('break')
        file_name = f'histos_{f_dmeson_species}_MBN{i}MCrec{fix_name}.json'
        if not exists(file_name):
            with open(file_name, 'w', encoding='utf-8') as f:
                json.dump(hist_list[i].tojson(), f, ensure_ascii=False, indent=4)

        if exists(file_name):
            with open(file_name, 'r') as f:
                json_data = json.loads(f.read())
                print(len(json_data))
                for item in json_data["arr"]:
                    if item['fName'] == 'hsDphiz':
                        sparse_dphiz.append(item)
                        print(item['fName'])
                        print(type(item))
                        print(len(item['fAxes']["arr"]))
                        # with open(f'test_{i}.json', 'w', encoding='utf-8') as t:
                        #     json.dump(item, t, ensure_ascii=False, indent=4)


djet_efficiency(mc_eff_file=run_mc_eff_file)
