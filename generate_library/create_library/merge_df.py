import pandas as pd


def merge_compound_feature_tables(compound_df, feature_df, mz_Da=0.002, mz_ppm=10):
    """
    Merge compound DataFrame with feature table based on m/z values.
    """

    # Create empty list to store matches
    matches = []

    # Iterate through each compound
    for _, compound in compound_df.iterrows():
        # Find all features within the mass tolerance
        matching_features = feature_df[
            (feature_df['m/z'] >= compound['t_mz'] - 0.02) &
            (feature_df['m/z'] <= compound['t_mz'] + 0.02)
            ]

        # If matches found, add them to the list
        if not matching_features.empty:
            for _, feature in matching_features.iterrows():
                match = compound.to_dict()
                match.update(feature.to_dict())
                match['mz_diff_ppm'] = abs(compound['t_mz'] - feature['m/z']) / compound['t_mz'] * 1e6
                match['mz_diff_Da'] = abs(compound['t_mz'] - feature['m/z'])
                match['selected'] = True
                match['discard_reason'] = ''
                match['ms2_explained_intensity'] = 0.0
                matches.append(match)

    # Create DataFrame from matches
    merged_df = pd.DataFrame(matches)

    # mz_diff
    merged_df = merged_df[(merged_df['mz_diff_ppm'] <= mz_ppm) | (merged_df['mz_diff_Da'] <= mz_Da)]
    merged_df = merged_df.reset_index(drop=True)

    return merged_df
