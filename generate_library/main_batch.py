import os

import pandas as pd

from feature_extraction import feature_extraction_single, plot_all_ms2
from cmpd import calculate_cmpd_mz
from filter_library import filter_library


def main_batch(file_dir, out_dir='.', data_collector='Minions'):
    """
    Process a batch of mzML files and csv files.
    """

    files = os.listdir(file_dir)
    files = [f for f in files if f.endswith('.mzML') and not f.startswith('.')]

    # create output directory
    os.makedirs(out_dir, exist_ok=True)

    metadata_dir = os.path.join(out_dir, 'metadata')
    os.makedirs(metadata_dir, exist_ok=True)

    all_library_df = pd.DataFrame()
    for f in files:
        print('Processing', f)
        mzml = os.path.join(file_dir, f)

        file_name = f.split('.')[0]
        csv = mzml.replace('.mzML', '.csv')

        # Extract features from mzML file
        feature_df = feature_extraction_single(file_path=mzml, save=False)

        # Calculate compound mz values
        print('Calculating compound mz values...')
        cmpd_df = calculate_cmpd_mz(csv)

        # Filter library
        print('Filtering library...')
        df, library_df = filter_library(cmpd_df, feature_df, data_collector, file_name)

        all_library_df = pd.concat([all_library_df, library_df])

        # Save the filtered library
        print('Saving the library...')
        out_df_path = os.path.join(metadata_dir, f'{file_name}_metadata.tsv')
        df.to_csv(out_df_path, sep='\t', index=False)

        # Plot all MS2 spectra
        print('Plotting all MS2 spectra...')
        plot_all_ms2(df, mzml, metadata_dir)

    all_library_df.to_csv(f'{out_dir}/all_library.tsv', sep='\t', index=False, na_rep='N/A')


if __name__ == '__main__':
    main_batch('/Users/shipei/Documents/projects/reverse_metabolomics/reverse_metabolomics/filter_library/input')

