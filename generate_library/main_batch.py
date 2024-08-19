import os

import pandas as pd

from feature_extraction import feature_extraction_single, plot_all_ms2, plot_all_eic
from cmpd import calculate_cmpd_mz
from filter_library import filter_library


def main_batch(file_dir,
               out_dir='.',
               data_collector='Minions',
               plot=True
               ):
    """
    Process a batch of mzML files and csv files.
    """

    files = os.listdir(file_dir)
    mzml_files = [f for f in files if f.endswith('.mzML') and not f.startswith('.')]
    csv_files = [f for f in files if f.endswith('.csv') and not f.startswith('.')]

    # create output directory
    os.makedirs(out_dir, exist_ok=True)
    # create metadata directory
    metadata_dir = os.path.join(out_dir, 'metadata')
    os.makedirs(metadata_dir, exist_ok=True)

    # Load the csv file
    print('Calculating compound mz values...')
    all_cmpd_df = pd.DataFrame()
    for _csv in csv_files:
        csv = os.path.join(file_dir, _csv)
        print('Loading', csv)
        cmpd_df = calculate_cmpd_mz(csv)
        all_cmpd_df = pd.concat([all_cmpd_df, cmpd_df])

    # Load the mzml file
    all_library_df = pd.DataFrame()
    for _mzml in mzml_files:
        print('Processing', _mzml)
        mzml = os.path.join(file_dir, _mzml)

        mzml_name = _mzml.split('.')[0]

        # Extract features from mzML file
        feature_df = feature_extraction_single(file_path=mzml, save=False)

        cmpd_df = all_cmpd_df[all_cmpd_df['unique_sample_id'] == _mzml].reset_index(drop=True).copy()

        # Filter library
        print('Filtering library...')
        df, library_df = filter_library(cmpd_df, feature_df, data_collector, mzml_name)

        all_library_df = pd.concat([all_library_df, library_df])

        # Save the filtered library
        print('Saving the library...')
        out_df_path = os.path.join(metadata_dir, f'{mzml_name}_metadata.tsv')
        df.to_csv(out_df_path, sep='\t', index=False)

        if plot:
            # Plot all MS2 spectra
            print('Plotting all MS2 spectra...')
            plot_all_ms2(df, mzml, metadata_dir)

            # Plot all EICs
            print('Plotting all EICs...')
            plot_all_eic(df, mzml, metadata_dir)

    all_library_df.to_csv(f'{out_dir}/all_library.tsv', sep='\t', index=False, na_rep='N/A')


if __name__ == '__main__':
    main_batch('/Users/shipei/Documents/projects/reverse_metabolomics/reverse_metabolomics/filter_library/input')

