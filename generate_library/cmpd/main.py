import pandas as pd
from .utils import neutralize_formula, calc_exact_mass


_adduct_pos = [
    {'name': '[M+H]+', 'm': 1, 'charge': 1, 'mass': 1.00727645223},
    {'name': '[M+Na]+', 'm': 1, 'charge': 1, 'mass': 22.989220702},
    {'name': '[M+K]+', 'm': 1, 'charge': 1, 'mass': 38.9631579064},
    {'name': '[M+NH4]+', 'm': 1, 'charge': 1, 'mass': 18.03382555335},
    # {'name': '[M-H2O+NH4]+', 'm': 1, 'charge': 1, 'mass': 0.023260869},
    {'name': '[M+H-H2O]+', 'm': 1, 'charge': 1, 'mass': -17.0032882318},
    {'name': '[M+H-2H2O]+', 'm': 1, 'charge': 1, 'mass': -35.01385291583},
    {'name': '[M+H-3H2O]+', 'm': 1, 'charge': 1, 'mass': -53.02441759986}
]

_adduct_neg = [
    {'name': '[M-H]-', 'm': 1, 'charge': 1, 'mass': -1.00727645223},
    {'name': '[M+Cl]-', 'm': 1, 'charge': 1, 'mass': 34.968304102},
    {'name': '[M+Br]-', 'm': 1, 'charge': 1, 'mass': 78.91778902},
    {'name': '[M+FA]-', 'm': 1, 'charge': 1, 'mass': 44.99710569137},
    {'name': '[M+Ac]-', 'm': 1, 'charge': 1, 'mass': 59.01275575583},
    {'name': '[M-H-H2O]-', 'm': 1, 'charge': 1, 'mass': -19.01784113626}
]


def calculate_cmpd_mz(cmpd_df_path, ion_mode='positive'):
    """
    Calculate the exact mass for each compound in the compound list
    """

    if ion_mode == 'positive':
        adduct_ls = _adduct_pos
    else:
        adduct_ls = _adduct_neg

    # load the compound list with USI and taxon filter
    cmpd_df = pd.read_csv(cmpd_df_path, low_memory=False)

    # neutralize the formula, deal with the charge (e.g., C5H5N+)
    cmpd_df['neutralized_formula'] = cmpd_df['formula'].apply(neutralize_formula)

    # calculate the exact mass
    cmpd_df['exact_mass'] = cmpd_df['neutralized_formula'].apply(calc_exact_mass)

    # Create a list to store the new rows
    new_rows = []

    # Iterate through each compound and adduct
    for _, compound in cmpd_df.iterrows():
        for adduct in adduct_ls:
            mz = compound['exact_mass'] + adduct['mass']
            new_row = compound.to_dict()
            new_row['t_mz'] = mz
            new_row['t_adduct'] = adduct['name']
            new_rows.append(new_row)

    # Create a new DataFrame from the list of new rows
    df = pd.DataFrame(new_rows)

    return df
