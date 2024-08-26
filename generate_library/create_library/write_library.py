import pandas as pd


def write_library(df, data_collector, file_name, ion_mode):
    """
    Write the filtered library to a file.
    """

    df = df[(df['selected']) & (~pd.isnull(df['best_MS2_scan_idx']))].reset_index(drop=True)

    charge = 1 if ion_mode == 'positive' else -1
    _ion_mode = 'Positive' if ion_mode == 'positive' else 'Negative'

    rows = []
    for _, row in df.iterrows():
        rows.append({
            'FILENAME': file_name + '.mzML',
            'SEQ': '*..*',
            'COMPOUND_NAME': row['name'],
            'MOLECULEMASS': row['t_mz'],
            'INSTRUMENT': 'Orbitrap',
            'IONSOURCE': 'LC-ESI',
            'EXTRACTSCAN': round(row['best_MS2_scan_idx']),
            'SMILES': row['SMILES'],
            'INCHI': row['inchi'],
            'INCHIAUX': row['inchi_adduct'],
            'CHARGE': charge,
            'IONMODE': _ion_mode,
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

    return pd.DataFrame(rows)
