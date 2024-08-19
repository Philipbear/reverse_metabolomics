import pandas as pd


def remove_smiles_with_empty_ms2(df):
    """
    Remove chemicals with empty MS2
    """
    unique_smiles = df['SMILES'].unique()
    for smiles in unique_smiles:
        mask = df['SMILES'] == smiles
        if pd.isnull(df.loc[mask, 'best_MS2_scan_idx']).all():
            df.loc[mask, 'selected'] = False
            df.loc[mask, 'discard_reason'] = 'No MS2 collected'

    return df


def remove_doubly_charged_ions(row):
    """
    Remove matches to doubly charged ions
    """
    if row['selected'] and row['charge'] == 2: #and len(row['isotopes']) > 2:
        row['selected'] = False
        row['discard_reason'] = 'Matched to doubly charged ion'
    return row


def remove_isotopes(row):
    """
    Remove matches to isotopes (the best MS2 scan index is empty)
    """
    if row['selected'] and row['is_isotope'] and row['best_MS2_scan_idx'] is None:
        row['selected'] = False
        row['discard_reason'] = 'Matched to isotope'
    return row
