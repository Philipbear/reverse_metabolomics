from .merge_df import merge_compound_feature_tables
from .filter_df import filter_df, filter_feature_df
from .write_library import write_library
from .write_mgf import write_mgf
from .summarize_df import summarize_df


def create_library(compound_df, feature_df, ion_mode, feature_intensity_threshold,
                   data_collector, file_name,
                   filter_library=True, ms2_explanation_cutoff=0.60,
                   metadata_dir=None,
                   write_individual_mgf=False):
    """
    Filter the library based on the compound and feature DataFrames.
    """

    # Filter the feature table
    feature_df = filter_feature_df(feature_df, feature_intensity_threshold)

    # Merge the compound and feature tables
    df = merge_compound_feature_tables(compound_df, feature_df)

    if filter_library:
        # Filter the merged dataframe
        df = filter_df(df, ion_mode, ms2_explanation_cutoff=ms2_explanation_cutoff)

    if write_individual_mgf:
        # Write individual MGF files
        write_mgf(df, file_name, metadata_dir)

    # summarize, number of peaks / isomers
    df = summarize_df(df)

    # Write the filtered library (dataframe to be uploaded to GNPS)
    library_df = write_library(df, data_collector, file_name, ion_mode)

    return df, library_df
