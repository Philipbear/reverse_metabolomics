from .merge_df import merge_compound_feature_tables
from .filter_df import filter_df
# from .write_library import write_library


def filter_library(compound_df, feature_df):
    """
    Filter the library based on the compound and feature DataFrames.
    """
    # Merge the compound and feature tables
    df = merge_compound_feature_tables(compound_df, feature_df)

    # Filter the merged DataFrame
    df = filter_df(df)

    # # Write the filtered library (dataframe to be uploaded to GNPS)
    # library_df = write_library(df)
    library_df = None

    return df, library_df
