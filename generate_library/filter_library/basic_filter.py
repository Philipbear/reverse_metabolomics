import pandas as pd
import numpy as np


def remove_smiles_with_empty_valid_ms2(df):
    """
    Remove chemicals with empty MS2
    """
    unique_smiles = df['SMILES'].unique()
    for smiles in unique_smiles:
        mask = df['SMILES'] == smiles

        _selected = df.loc[mask, 'selected']
        _best_MS2_scan_idx = df.loc[mask, 'best_MS2_scan_idx']
        _bool = _selected & pd.notnull(_best_MS2_scan_idx)
        if not _bool.any():
            df.loc[mask, 'selected'] = False
            df.loc[mask, 'discard_reason'] = np.where(
                df.loc[mask, 'discard_reason'] == '', 'No valid MS2',
                df.loc[mask, 'discard_reason'] + '; No valid MS2'
            )

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
