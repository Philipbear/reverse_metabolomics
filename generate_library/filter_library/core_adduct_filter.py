import pandas as pd
import numpy as np


def core_adduct_filter(df, core_adduct_ls=None, rt_tol=0.025):
    """
    For each SMILES and grouped RT, check if core adducts are present
    """
    if core_adduct_ls is None:
        core_adduct_ls = ['[M+H]+', '[M+NH4]+', '[M+H-H2O]+', '[M+H-2H2O]+', '[M+H-3H2O]+']

    # Find all unique SMILES
    unique_smiles = df['SMILES'].unique()

    for smiles in unique_smiles:
        # Get all rows for this SMILES
        smiles_df = df[df['SMILES'] == smiles]

        # Group RTs for this SMILES
        rt_groups = group_rts(smiles_df['RT'].unique(), rt_tol)

        for rt_group in rt_groups:
            # Get all rows for this SMILES and RT group
            rt_mask = smiles_df['RT'].isin(rt_group)
            rt_df = smiles_df[rt_mask]

            # Check if any core adduct is present
            core_adduct_present = rt_df['t_adduct'].isin(core_adduct_ls).any()

            if not core_adduct_present:
                # If no core adduct is present, update the original dataframe
                update_mask = (df['SMILES'] == smiles) & df['RT'].isin(rt_group)
                df.loc[update_mask, 'selected'] = False

                # Update discard reason
                df.loc[update_mask, 'discard_reason'] = np.where(
                    df.loc[update_mask, 'discard_reason'] == '',
                    'No core adduct',
                    df.loc[update_mask, 'discard_reason'] + ';No core adduct'
                )

    return df


def group_rts(rts, rt_tol):
    """Group RTs that are within the tolerance of each other"""
    sorted_rts = sorted(rts)
    groups = []
    current_group = [sorted_rts[0]]

    for rt in sorted_rts[1:]:
        if rt - current_group[0] <= rt_tol:
            current_group.append(rt)
        else:
            groups.append(current_group)
            current_group = [rt]

    groups.append(current_group)
    return groups
