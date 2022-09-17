"""
/**************************************************************************
 * Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
"""

# -------------------------------------------------------------------------
# Author: Auro Mohanty
# Utrecht University
# auro.mohanty@cern.ch 
# -------------------------------------------------------------------------
import uproot

from run_settings import run_is_eff


def signal_extraction_sb(
        data_path: str = "",
        is_eff: bool = run_is_eff,
        eff_file: str = "",  # file after having been processed by djet_efficiency.py
        is_ref: bool = 0,
        ref_file: str = "",
        post_fix: bool = 0,
        list_name: str = "",  # Cut
        out: str = "signal_extraction",
        save: bool = 1,
        is_more_files: bool = 0,
        prod: str = "kl",
        is_prefix: str = 0,
        #
        # raw yield systematics: multi-trial
        bound_sigma: int = 0,
        f_sigma_factor: float = 1,
        fixed_mass: bool = 0,
):
    # efficiency
    if is_eff:
        file_eff = uproot.open(eff_file)
        print(file_eff)
        for key in file_eff.keys():
            print(key)
        # h_eff = file_eff["hEff_reb"]
        # print(h_eff)


signal_extraction_sb(is_eff=run_is_eff)
