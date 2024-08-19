import os
from feature_extraction import feature_extraction_single
from cmpd import calculate_cmpd_mz
from filter_library import filter_library


def main_single(mzml, out_dir='.',
                data_collector='Minions'):

    file_name = mzml.split('/')[-1].split('.')[0]
    csv = mzml.replace('.mzML', '.csv')

    # Extract features from mzML file
    feature_df = feature_extraction_single(file_path=mzml, save=False)

    # Calculate compound mz values
    print('Calculating compound mz values...')
    cmpd_df = calculate_cmpd_mz(csv)

    # Filter library
    print('Filtering library...')
    df, library_df = filter_library(cmpd_df, feature_df, data_collector, file_name)

    # Save the filtered library
    print('Saving the library...')
    os.makedirs(out_dir, exist_ok=True)
    df.to_csv(f'{out_dir}/{file_name}_metadata.tsv', sep='\t', index=False)
    library_df.to_csv(f'{out_dir}/{file_name}.tsv', sep='\t', index=False, na_rep='N/A')


if __name__ == '__main__':
    main_single('/Users/shipei/Documents/projects/reverse_metabolomics/reverse_metabolomics/filter_library/input/AP_68.mzML')

