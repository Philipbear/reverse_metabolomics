import pandas as pd
import os


def write_mgf(df, file_name, out_dir):
    """
    Write the filtered library to a file.
    """

    df = df[~pd.isnull(df['best_MS2_scan_idx'])].reset_index(drop=True)

    out_path = os.path.join(out_dir, f'{file_name}_allMS2.mgf')

    with open(out_path, 'w') as f:
        for _, row in df.iterrows():
            f.write('BEGIN IONS\n')
            f.write(f'SELECTED={row["selected"]}\n')
            f.write(f'COMPOUND={row["compound_name"]}\n')
            f.write(f'PEPMASS={row["t_mz"]}\n')
            f.write(f'ION={row["t_adduct"]}\n')
            f.write(f'FILENAME={file_name}.mzML\n')
            f.write(f'SCANS={round(row["best_MS2_scan_idx"])}\n')

            mzs = row['MS2'][:, 0]
            intensities = row['MS2'][:, 1]
            for mz, intensity in zip(mzs, intensities):
                f.write(f'{mz} {intensity}\n')
            f.write('END IONS\n')

    return
