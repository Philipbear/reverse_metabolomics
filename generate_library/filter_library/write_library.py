import pandas as pd


def write_library(df, data_collector, file_name):
    """
    Write the filtered library to a file.
    """

    df = df[(df['selected']) & (pd.isnull(df['all_MS2_scan_idx']) == False)].reset_index(drop=True)

    rows = []
    for _, row in df.iterrows():
        rows.append({
            'FILENAME': file_name + '.mzML',
            'SEQ': '*..*',
            'COMPOUND_NAME': row['compound_name'],
            'MOLECULEMASS': row['t_mz'],
            'INSTRUMENT': 'Orbitrap',
            'IONSOURCE': 'LC-ESI',
            'EXTRACTSCAN': round(row['best_MS2_scan_idx']),
            'SMILES': row['SMILES'],
            'INCHI': None,
            'INCHIAUX': None,
            'CHARGE': 1,
            'IONMODE': 'Positive',
            'PUBMED': None,
            'ACQUISITION': 'Crude',
            'EXACTMASS': row['exact_mass'],
            'DATACOLLECTOR': data_collector,
            'ADDUCT': row['t_adduct'],
            'CASNUMBER': None,
            'PI': 'Pieter Dorrestein',
            'LIBQUALITY': 1,
            'GENUS': None,
            'SPECIES': None,
            'INTEREST': None,
            'STRAIN': None
        })


    library_df = pd.DataFrame(rows)

    return library_df
