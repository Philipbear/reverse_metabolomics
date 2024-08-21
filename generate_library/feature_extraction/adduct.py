import numpy as np


def peak_peak_correlation(roi1, roi2):
    """
    A function to find the peak-peak correlation between two rois.

    roi: ROI object
    """

    # find the common scans in the two rois
    common_scans = np.intersect1d(roi1.scan_idx_seq, roi2.scan_idx_seq)

    if len(common_scans) < 4:
        return 0.0

    # find the intensities of the common scans in the two rois
    int1 = roi1.int_seq[np.isin(roi1.scan_idx_seq, common_scans)]
    int2 = roi2.int_seq[np.isin(roi2.scan_idx_seq, common_scans)]

    # calculate the correlation
    # if all values are same, return 1
    if np.all(int1 == int1[0]) or np.all(int2 == int2[0]):
        return 1.0

    pp_cor = np.corrcoef(int1, int2)[0, 1]
    return pp_cor


def annotate_adduct(d):
    """
    A function to annotate adduct
    """
    mz_tol = d.params.mz_tol_ms1
    rt_tol = d.params.peak_cor_rt_tol

    # sort rois by mz
    d.rois.sort(key=lambda x: x.mz)
    roi_to_label = np.ones(len(d.rois), dtype=bool)

    # isotopes, in-source fragments
    for idx, r in enumerate(d.rois):
        if r.is_isotope or r.is_in_source_fragment:
            roi_to_label[idx] = False

    if d.params.ion_mode.lower() == "positive":
        default_single_adduct = "[M+H]+"
        default_double_adduct = "[M+2H]2+"
    elif d.params.ion_mode.lower() == "negative":
        default_single_adduct = "[M-H]-"
        default_double_adduct = "[M-2H]2-"

    adduct_list = _pop_adduct_list(_adduct_pos, _adduct_neg, d.params.ion_mode.lower())

    # find adducts by assuming the current roi is [M+H]+ ([M-H]- in negative mode)
    group_no = 1
    for idx, r in enumerate(d.rois):

        if not roi_to_label[idx]:
            continue

        if r.charge_state != 1:
            continue

        mol_mass = _calc_exact_mol_mass(r.mz, d.params.ion_mode.lower())

        # for every possible adduct
        for i, adduct in enumerate(adduct_list):
            m = (mol_mass * adduct['m'] + adduct['mass']) / adduct['charge']
            _v = np.logical_and(roi_to_label, np.abs(d.roi_rt_seq - r.rt) < rt_tol)

            v = np.where(np.logical_and(_v, np.abs(d.roi_mz_seq - m) < mz_tol))[0]

            if len(v) == 0:
                continue

            if len(v) > 1:
                # select the one with the lowest RT difference
                v = v[np.argmin(np.abs(d.roi_rt_seq[v] - r.rt))]
            else:
                v = v[0]

            if peak_peak_correlation(r, d.rois[v]) > d.params.ppr:
                roi_to_label[v] = False
                d.rois[v].adduct_type = adduct['name']
                d.rois[v].adduct_parent_roi_id = r.id
                d.rois[v].adduct_group_no = group_no
                r.adduct_child_roi_id.append(d.rois[v].id)

        if len(r.adduct_child_roi_id) > 0:
            roi_to_label[idx] = False
            r.adduct_type = default_single_adduct
            r.adduct_group_no = group_no
            group_no += 1

    for r in d.rois:
        if r.adduct_type is None:
            r.adduct_type = default_single_adduct if r.charge_state == 1 else default_double_adduct


def _calc_exact_mol_mass(mz, ion_mode):
    if ion_mode == 'positive':
        return mz - 1.00727645223
    else:
        return mz + 1.00727645223


def _pop_adduct_list(_adduct_pos, _adduct_neg, ion_mode):
    if ion_mode == 'positive':
        return _adduct_pos[1:]
    else:
        return _adduct_neg[1:]


_adduct_pos = [
    {'name': '[M+H]+', 'm': 1, 'charge': 1, 'mass': 1.00727645223},
    {'name': '[M+Na]+', 'm': 1, 'charge': 1, 'mass': 22.989220702},
    {'name': '[M+K]+', 'm': 1, 'charge': 1, 'mass': 38.9631579064},
    {'name': '[M+NH4]+', 'm': 1, 'charge': 1, 'mass': 18.03382555335},
    {'name': '[M+H-H2O]+', 'm': 1, 'charge': 1, 'mass': -17.0032882318},
    {'name': '[M+H-2H2O]+', 'm': 1, 'charge': 1, 'mass': -35.01385291583},
    {'name': '[M+H-3H2O]+', 'm': 1, 'charge': 1, 'mass': -53.02441759986}
]

_adduct_neg = [
    {'name': '[M-H]-', 'm': 1, 'charge': 1, 'mass': -1.00727645223},
    {'name': '[M+Cl]-', 'm': 1, 'charge': 1, 'mass': 34.968304102},
    {'name': '[M+Br]-', 'm': 1, 'charge': 1, 'mass': 78.91778902},
    {'name': '[M+FA]-', 'm': 1, 'charge': 1, 'mass': 44.99710569137},
    {'name': '[M+Ac]-', 'm': 1, 'charge': 1, 'mass': 59.01275575583},
    {'name': '[M-H-H2O]-', 'm': 1, 'charge': 1, 'mass': -19.01784113626}
]
