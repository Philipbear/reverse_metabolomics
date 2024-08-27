from msbuddy import assign_subformula


def filter_by_ms2_explanation(row, explanation_cutoff=0.60):
    """
    Calculate the MS2 explanation
    """
    if row['selected'] and row['MS2'] is not None:
        explained_intensity = 0.0
        subformla_list = assign_subformula(row['MS2'][:, 0],
                                           precursor_formula=row['neutralized_formula'],
                                           adduct=row['t_adduct'],
                                           ms2_tol=0.01, ppm=False)
        for subformula in subformla_list:
            explained_intensity += row['MS2'][:, 1][subformula.idx] if subformula.subform_list else 0.0
        row['ms2_explained_intensity'] = explained_intensity / row['MS2'][:, 1].sum()

        if row['ms2_explained_intensity'] < explanation_cutoff:
            row['selected'] = False
            row['discard_reason'] = 'MS2 explanation below cutoff'

    return row
