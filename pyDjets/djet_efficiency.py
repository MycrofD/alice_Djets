import uproot
import numpy as np
from central_config import f_dmeson_species, NDMC
from run_settings import run_eff_file, run_is_prefix, run_is_postfix


def djet_efficiency(
    eff_file: str = run_eff_file,
    is_prefix: bool = run_is_prefix,
    is_postfix: bool = run_is_postfix,
    fix_name: str = "",
):
    root_file = uproot.open(eff_file)
    directory_file = root_file["DmesonsForJetCorrelations"]

    hist_list = []
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

        print(vars(hist_list[i]))
        print(hist_list[i].members)
    # hist_list[0] is [hstat, hCentDjet, ...]
    print("break 0=======")
    shs_dphiz = np.array(hist_list[0].bases)
    for item in shs_dphiz[0]:
        print("====-----------------")
        print(vars(item))
        # print(item.name)
    print("break 1=======")
    t = shs_dphiz[0][-1].bases
    print(type(t))
    print(len(t))
    for item in t:
        print(vars(item))
        print("break --- =======")
        print(vars(item.bases[0]))
        print(item.bases[0])
        print("break --- =======")
        print(item.bases[0].members)
    print("break 2=======")


djet_efficiency(eff_file=run_eff_file)
