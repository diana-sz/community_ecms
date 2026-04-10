#!/usr/bin/env python
# coding: utf-8


import cobra
#import efmtool
import numpy as np
#import pandas as pd

cobra_config = cobra.Configuration()
cobra_config.solver = "cplex"
cobra_config.tolerance = 1e-9


def flip_reverse_reactions(model, lb_threshold=-999.0, ub_threshold=0.0, verbose=False):
    """
    Flip reactions in a COBRApy model that only allow reverse flux
    (e.g. lower_bound ~ -1000 and upper_bound = 0) so they become forward reactions.

    Parameters
    ----------
    model : cobra.Model
        The COBRApy model
    lb_threshold : float
        Threshold to detect "very negative" lower bounds (default: -999)
    ub_threshold : float
        Threshold for upper bound (default: 0)
    verbose : bool
        Print flipped reactions

    Returns
    -------
    flipped_rxns : list
        List of reaction IDs that were flipped
    """

    flipped_rxns = []

    for rxn in model.reactions:
        lb = rxn.lower_bound
        ub = rxn.upper_bound

        # Detect reverse-only reactions
        if lb <= lb_threshold and ub == ub_threshold:
            # Flip stoichiometry
            new_stoich = {met: -coeff for met, coeff in rxn.metabolites.items()}
            rxn.add_metabolites(new_stoich, combine=False)

            # Set bounds to forward direction
            rxn.lower_bound = 0.0
            rxn.upper_bound = -lb  # e.g. 1000

            flipped_rxns.append(rxn.id)

            if verbose:
                print(f"Flipped reaction: {rxn.id}")

    return flipped_rxns





bm_rxns = {"Geobacter": "agg_GS13m", "Rhodoferax": "BIO_Rfer3"}
c_sources = ["EX_ac_e", "EX_mal-L_e", "EX_cit_e", "EX_fum_e"]


dimes = {"Rhodoferax": {
            "mEmC_13": {"uptakes": ["EX_nh4_e", "EX_fe2_e", "EX_fum_e",
                                  "EX_so4_e", "EX_pi_e", "EX_o2_e"],
                "secretion": ["EX_mal-L_e",  "EX_co2_e"]},
            "mEmC_24": {"uptakes": ["EX_nh4_e", "EX_fe2_e", "EX_cit_e",
                                  "EX_so4_e", "EX_pi_e", "EX_o2_e"],
                "secretion": ["EX_mal-L_e",  "EX_co2_e"]},
            "mEmC_5": {"uptakes": ["EX_nh4_e", "EX_fe2_e", "EX_mal-L_e",
                                  "EX_so4_e", "EX_pi_e", "EX_o2_e"],
                "secretion": ["EX_co2_e"]},
            "mEmC_6": {"uptakes": ["EX_nh4_e", "EX_fe2_e", "EX_mal-L_e", 
                                  "EX_so4_e", "EX_pi_e", "EX_o2_e",
                                  "EX_co2_e"],
                "secretion": ["EX_orot_e"]},
        "mE_1": {"uptakes": ["EX_nh4_e", "EX_ac_e", "EX_fe3_e",
                            "EX_so4_e", "EX_pi_e", "EX_o2_e"],
                "secretion": ["EX_co2_e"]},
},
        "Geobacter": {
            "mEmC_12": {"uptakes": ["EX_nh4_e", "EX_fe3_e", "EX_so4_e", "EX_pi_e",
                                  "EX_mal-L_e"],
                      "secretion": ["EX_fe2_e", "EX_co2_e", "EX_h2s_e"]},
            "mEmC_34": {"uptakes": ["EX_n2_e", "EX_fe3_e", "EX_so4_e", "EX_pi_e",
                                  "EX_mal-L_e"],
                      "secretion": ["EX_fe2_e", "EX_co2_e", "EX_h2s_e"]},
            "mEmC_56": {"uptakes": ["EX_n2_e", "EX_fe3_e", "EX_so4_e",
                                  "EX_pi_e", "EX_cit_e"],
                      "secretion": ["EX_fe2_e", "EX_mal-L_e", "EX_co2_e",
                                    "EX_h2s_e"]},
            "mE_1": {"uptakes": ["EX_n2_e", "EX_fe3_e", "EX_so4_e",
                              "EX_pi_e", "EX_mal-L_e"],
                  "secretion": ["EX_fe2_e", "EX_ac_e", "EX_co2_e",
                                "EX_h2s_e"]},
    }}

# from their github
def open_exchanges(model, organism):
    
    if organism == 'Geobacter':
        allowed = [#'h_e',
             # 'fe2_e',
             # 'co2_e',
             'so4_e',
             'pi_e',
             'mg2_e',
             'k_e',
             #'h2o_e',
             #'fe3_e',
             #'zn2_e',
             #'ss_e',
             #'s_e',
             #'ni2_e',
             #'na1_e',
             #'n2_e',           
             #'mobd_e',
             #'mn2_e',
             # 'h2s_e',
             # 'h2_e',
             #'cu2_e',
             #'cobalt2_e',
             #'cl_e',
             #'cd2_e',
             'ca2_e']
    elif organism == 'Rhodoferax':
        allowed = [
            # '2ddglcn_e',
             # 'ac_e',
             # 'arab-L_e',
             # 'buts_e',
             # 'bz_e',
             #'ca2_e',
             #'cd2_e',
             # 'cellb_e',
             # 'cit_e',
             # 'co2_e',
             #'cobalt2_e',
             #'cro4_e',
             #'cu2_e',
             # 'eths_e',
             # 'etoh_e',
             # 'fe2_e',
             'fe3_e',
             # 'fru_e',
             # 'fum_e',
             # 'glc_e',
             # 'glyclt_e',
             'h_e',
             'h2_e',
             'h2o_e',
             # 'hexs_e',
             # 'istnt_e',
             #'k_e',
             # 'lac-L_e',
             # 'mal-L_e',
             #'mg2_e',
             #'mn2_e',
             #'mobd_e',
             #'na1_e',
             'nh4_e',
             #'ni2_e',
             # 'no2_e',
             # 'no3_e',
             # 'o2_e',
             'pi_e',
             # 'ppa_e',
             # # 'ppn_e',
             # 'pyr_e',
             # 'rib-D_e',
             'so4_e',
             # 'succ_e',
             # 'sula_e',
             #'tsul_e',
             #'zn2_e',
            ]
    elif organism == "Shewanella":
        allowed = ['ac_e',
                # 'akg_e',
                'arsna_e',
                'arsni2_e',
                'ca2_e',
                'cl_e',
                # 'co2_e',
                'cobalt2_e',
                'cobalt3_e',
                'CrOH3_e',
                'cro4_e',
                'cu2_e',
                # 'dms_e',
                # 'dmso_e',
                # 'etoh_e',
                # 'fe2_e',
                'fe3_e',
                # 'for_e',
                # 'fum_e',
                # 'glyclt_e',
                # 'glyc-R_e',
                'h_e',
                'h2_e',
                'h2o_e',
                'h2o2_e',
                'h2s_e',
                # 'hdca_e',
                'hg2_e',
                # 'inoshp_e',
                # 'inospp1_e',
                'k_e',
                # # 'lac-D_e',
                # 'lac-L_e',
                # 'mal-L_e',
                'mg2_e',
                'mn2_e',
                'mn4o_e',
                'mobd_e',
                'na1_e',
                'nh4_e',
                'ni2_e',
                # 'no2_e',
                # 'no3_e',
                # 'o2_e',
                # 'ocdca_e',
                'pi_e',
                # 'ppa_e',
                # 'pyr_e',
                'so3_e',
                'so4_e',
                # 'succ_e',
                'tsul_e',
                # 'ttdca_e',
                'tttnt_e',
                'wo4_e',]

    # we close demand reactions here, as the other functions do not care about the demands
    # close the demand reactions in Geobacter

    if organism == "Geobacter":
        to_remove = [rxn.id for rxn in model.reactions if rxn.id.startswith("DM_")]

        model.remove_reactions(to_remove)

    return allowed


def implement_dime(model, allowed, dime, allow_secretions=False):
    # first, allow all exchanges as in the paper and turn off all secretions
    mu = model.slim_optimize() 
    print(f" * Growth rate before implementing DiME {mu:.2f}")

    to_remove = []
    for ex in model.exchanges:
        if ex.id in [ "EX_h2o_e",  "EX_h2_e",  "EX_h_e" ]:
            ex.bounds = (-1000, 1000)
            continue

        # dime should overwrite the allowed
        if ex.id in dime["secretion"]:
            ex.bounds = (0, 1000)
        elif (ex.id.replace("EX_", "") in allowed) or (ex.id in dime["uptakes"]):
            ex.bounds = (-1000, 0)
        else:
            ex.bounds = (0, 1000)
            if not allow_secretions:
                to_remove.append(ex)

    model.remove_reactions(to_remove)
    flip_reverse_reactions(model)
    cobra.manipulation.delete.prune_unused_metabolites(model)

    mu = model.slim_optimize() 
    print(f" * Growth rate after implementing DiME {mu:.2f}")

    return mu



def find_inactive_rxns_fva(model, min_obj=1, TOL=1e-10):
    fva = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=min_obj)
    off = []
    upper_zero = []
    lower_zero = []
    rxn_dict = {r.id: r for r in model.reactions}

    for rxn_id, row in fva.iterrows():
        r = rxn_dict[rxn_id]
        min_flux = row["minimum"]
        max_flux = row["maximum"]

        if (abs(min_flux) < TOL) and (abs(max_flux) < TOL):
            off.append(rxn_id)
            continue

        if min_flux > 0:#-TOL:
            lower_zero.append(rxn_id)
        if max_flux < TOL:
            upper_zero.append(rxn_id)

    return(off, lower_zero, upper_zero)



reduced_models = {"Rhodoferax": {}, "Geobacter": {}}
for allow_secretions in [True, False]:
    for modelname in reduced_models:
        print(f"\n***** {modelname} *****")
        for dime_number, dime in dimes[modelname].items():
            print(f"Testing DiME {dime_number}")
    
            model = cobra.io.load_matlab_model(f"models/{modelname}.mat")
            allowed = open_exchanges(model, modelname)
            flip_reverse_reactions(model)
    
            to_remove = []
            for rxn in model.reactions:
                if abs(rxn.lower_bound) < 1e-15 and abs(rxn.upper_bound) < 1e-15:
                    print(f"{rxn.id} removed\n")
                    to_remove.append(rxn)
            model.remove_reactions(to_remove)
            cobra.manipulation.delete.prune_unused_metabolites(model)
    
            n_reactions = len(model.reactions)
        
            mu = model.slim_optimize()
            if not dime == "default":
                mu = implement_dime(model, allowed, dime, allow_secretions)
                suffix = ["nex", "ex"][allow_secretions]
                cobra.io.write_sbml_model(model, f"models/{modelname}_{dime_number}_{suffix}.xml")

            if np.isnan(mu):
                print(f" * no feasible growth with dime {dime_number}")
                continue
    
            rxn_dict = {r.id: r for r in model.reactions}

            # n = 0
            # min_obj = 1
            # while True:
            #     off, lower_zero, upper_zero = find_inactive_rxns_fva(model, min_obj=min_obj)
            #     if len(off) == 0 or n > 1:
            #         break
    
            #     n += 1
    
            #     model.remove_reactions([rxn_dict[rid] for rid in off])
            #     print(f" * FVA{n}: {len(off)} reaction(s) removed")

            #     for rxn_id in lower_zero:
            #         r = rxn_dict[rxn_id]
            #         r.lower_bound = 0

            #     for rxn_id in upper_zero:
            #         r = rxn_dict[rxn_id]
            #         r.upper_bound = 0
    
            # print(f"  ** Growth rate with FVA{n} constraints {model.slim_optimize():.2f}")       
            # print(f"  ** {len(model.reactions)}/{n_reactions} reactions in final model")

            # flip_reverse_reactions(model)
            # cobra.manipulation.delete.prune_unused_metabolites(model)
            
            solution = model.optimize()
            for c_source in c_sources:
                if c_source in model.medium:
                    bm_yield = solution.fluxes[bm_rxns[modelname]]/solution.fluxes[c_source]
                    print(f"Biomass yield: {bm_yield:.3f}")
                    break

            # reduced_models[modelname][dime_number] = model
            # cobra.io.write_sbml_model(model, f"{modelname}_{dime_number}_reduced_{suffix}_{min_obj}.xml")

