import pandas as pd
from feature_extraction import feature_extraction_single
from cmpd import calculate_cmpd_mz
from filter_library import filter_library


def main_single(mzml, csv):

    # Extract features from mzML file
    feature_df = feature_extraction_single(file_path=mzml, save=False)

    # Calculate compound mz values
    print('Calculating compound mz values...')
    cmpd_df = calculate_cmpd_mz(csv)

    # Filter library
    print('Filtering library...')
    df, library_df = filter_library(cmpd_df, feature_df)

    df.to_csv('filtered_library.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main_single('/Users/shipei/Documents/projects/reverse_metabolomics/reverse_metabolomics/filter_library/input/AP_68.mzML',
                '/Users/shipei/Documents/projects/reverse_metabolomics/reverse_metabolomics/filter_library/input/AP_68.csv')

