import pandas as pd
from .basic_filter import remove_smiles_with_empty_ms2, remove_doubly_charged_ions, remove_isotopes
from .ms2_explanation_filter import calc_ms2_explanation


def filter_df(df):
    """
    Filter the merged DataFrame
    """

    # remove smiles with empty MS2 at all
    df = remove_smiles_with_empty_ms2(df)

    # if matched to a doubly charged ion, remove the match
    df = df.apply(remove_doubly_charged_ions, axis=1)

    # if matched to an isotope, remove the match
    df = df.apply(remove_isotopes, axis=1)

    # calculate the MS2 explained intensity
    print('Calculating MS2 explained intensity...')
    df = df.apply(calc_ms2_explanation, axis=1)


    return df




