"""
This script is used to filter the library generated by MZmine3 workflow.
Filter the library based on core adduct list, ion dependency, modcos match, and isobaric mass check
@Author: Shipei Xing, Ipsita Mohanty
@Date: May 1, 2024
"""

import os

import numpy as np
import pandas as pd
from matchms import Spectrum
from matchms.similarity import ModifiedCosine, CosineGreedy
from molmass import Formula
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def generate_library_df(library_mgf, name_sep='_'):
    """
    Generate metadata dataframe for the mgf file
    name_sep: separator for the compound name. eg, 'Phe_CA' -> '_' is the separator
    """
    mgf_file = os.path.join('input', library_mgf)
    with open(mgf_file, 'r') as file:
        spectrum_list = []
        db_idx = 1
        for line in file:
            # empty line
            _line = line.strip()  # remove leading and trailing whitespace
            if not _line:
                continue
            elif line.startswith('BEGIN IONS'):
                spectrum = {}
                # initialize spectrum
                mz_list = []
                intensity_list = []
            elif line.startswith('END IONS'):
                if len(mz_list) == 0:
                    continue
                spectrum['mz_ls'] = mz_list
                spectrum['intensity_ls'] = intensity_list
                spectrum['db_idx'] = db_idx
                db_idx += 1
                spectrum_list.append(spectrum)
                continue
            else:
                # if line contains '=', it is a key-value pair
                if '=' in _line:
                    # split by first '='
                    key, value = _line.split('=', 1)
                    spectrum[key] = value
                else:
                    # if no '=', it is a spectrum pair
                    this_mz, this_int = _line.split()
                    try:
                        mz_list.append(float(this_mz))
                        intensity_list.append(float(this_int))
                    except:
                        continue

    df = pd.DataFrame(spectrum_list)

    # split adduct by '[' and ']', get the middle part
    df['_ADDUCT'] = df['ADDUCT'].apply(lambda x: x.split('[')[1].split(']')[0] if '[' in x else x)

    # split name by name_sep
    df['_NAME'] = df['NAME'].apply(lambda x: x.split(' (Chimeric')[0] if ' (Chimeric' in x else x)
    df['_NAME'] = df['_NAME'].apply(lambda x: x.split('_NCE')[0] if '_NCE' in x else x)
    df['NAME_1'] = df['_NAME'].apply(lambda x: x.split(name_sep, 1)[0] if name_sep in x else None)
    df['NAME_2'] = df['_NAME'].apply(lambda x: x.split(name_sep, 1)[1] if name_sep in x else None)
    df['conjugate'] = df['NAME_2'].apply(lambda x: True if x is not None else False)

    return df


def select_library(df, cmpd_df_dict, df_base_name, core_adduct_ls=None, rt_tol=2,
                   ms2_tol_da=0.02, prec_intensity_cutoff=10,
                   modcos_score_cutoff=0.6, modcos_peak_cutoff=4, cos_score_cutoff=0.95,
                   write_df=False):
    """
    Filter the library based on core adduct list, ion dependency, modcos match, and isobaric mass check
    if observing one precursor existence beyond prec_intensity_cutoff, remove ones with prec intensity 0.
    modcos: either pass the score or the number of matched peaks
    """
    df['selected'] = [True] * df.shape[0]
    df['discard_reason'] = [None] * df.shape[0]

    # convert to float
    df['RTINSECONDS'] = df['RTINSECONDS'].astype(float)
    df['PEPMASS'] = df['PEPMASS'].astype(float)
    # bin the RT
    df['_RT'] = pd.cut(df['RTINSECONDS'], bins=range(0, int(df['RTINSECONDS'].max()) + 1, rt_tol))
    # bin the PEPMASS
    df['_PEPMASS'] = df['PEPMASS'].apply(lambda x: round(x, 3))

    ##########################################
    # filter spectra, indicated by precursor existence
    df['cmpd_1_prec_int'] = [None] * df.shape[0]
    df['cmpd_2_prec_int'] = [None] * df.shape[0]
    df['max_prec_int'] = [0] * df.shape[0]
    discarded_scan_ls = precursor_check(df, cmpd_df_dict,
                                        ms2_tol_da=ms2_tol_da,
                                        prec_intensity_cutoff=prec_intensity_cutoff)
    df.loc[df['db_idx'].isin(discarded_scan_ls), 'selected'] = False
    df.loc[df['db_idx'].isin(discarded_scan_ls), 'discard_reason'] = 'precursor_check'

    print(f'{df.shape[0]} spectra in the library')
    print('After precursor existence check:')
    print(f'{df["selected"].sum()} spectra remaining')
    print(f'{len(df["NAME"][df["selected"]].unique())} compounds')

    ##########################################
    # core adduct check
    discarded_scan_ls, cmpd_adduct_to_be_removed_dict = core_adduct_check(df, core_adduct_ls=core_adduct_ls,
                                                                          ms2_tol_da=ms2_tol_da,
                                                                          modcos_score_cutoff=modcos_score_cutoff,
                                                                          modcos_peak_cutoff=modcos_peak_cutoff)
    df.loc[df['db_idx'].isin(discarded_scan_ls), 'selected'] = False
    df.loc[df['db_idx'].isin(discarded_scan_ls), 'discard_reason'] = 'core_adduct_check'

    print('After adduct check, ion dependency, mod cos match:')
    print(f'{df["selected"].sum()} spectra remaining')
    print(f'{len(df["NAME"][df["selected"]].unique())} compounds')

    ##########################################
    # remove spectra, indicated by cmpd_adduct_to_be_removed_dict
    df = isobaric_mass_adduct_check(df, cmpd_adduct_to_be_removed_dict)

    print('After isobaric mass filter (unique adduct existence):')
    print(f'{df["selected"].sum()} spectra selected')
    print(f'{len(df["NAME"][df["selected"]].unique())} compounds')

    ##########################################
    # remove almost identical spectra
    df = remove_identical_spectra(df, ms2_tol_da=ms2_tol_da, cos_score_cutoff=cos_score_cutoff)

    print('After removing almost identical spectra:')
    print(f'{df["selected"].sum()} spectra selected')
    print(f'{len(df["NAME"][df["selected"]].unique())} compounds')

    ##########################################
    # remove doubly charged adduct
    _doubly_charged_adducts = df['ADDUCT'].str.contains(r'\][\s]*[+]?2|2[+]', regex=True)
    df.loc[_doubly_charged_adducts, 'selected'] = False
    df.loc[_doubly_charged_adducts, 'discard_reason'] = 'doubly_charged_adduct'

    print('After removing doubly charged adducts:')
    print(f'{df["selected"].sum()} spectra remaining')
    print(f'{len(df["NAME"][df["selected"]].unique())} compounds')

    if write_df:
        file_path = os.path.join('output', 'library_metadata')
        file_path = os.path.join(file_path, f'{df_base_name}_metadata.tsv')
        df.to_csv(file_path, sep='\t', index=False)

    # remove cols
    df = df.drop(['_ADDUCT', '_NAME', 'db_idx', '_RT', 'NAME_1', 'NAME_2', 'conjugate', '_PEPMASS', 'cmpd_1_prec_int',
                  'cmpd_2_prec_int', 'max_prec_int', 'discard_reason'], axis=1)

    return df


def remove_identical_spectra(df, ms2_tol_da=0.02, cos_score_cutoff=0.95):
    """
    Remove almost identical spectra
    """
    # sort df by NAME and ADDUCT
    df = df.sort_values(['NAME', 'ADDUCT'])

    # Create a copy of the dataframe to store updates
    df_updated = df.copy()

    # Group by NAME and ADDUCT
    for (name, adduct), group in df.groupby(['NAME', 'ADDUCT']):
        # Sort group by precursor_purity from high to low and reset index
        group = group.sort_values('PRECURSOR_PURITY', ascending=False).reset_index()

        # Iterate through the group
        for i in range(1, len(group)):
            if group.loc[i, 'selected']:
                for j in range(i):
                    if group.loc[j, 'selected']:
                        # Compare spectra
                        score = compare_spectra_cos(group.iloc[j], group.iloc[i], ms2_tol_da)
                        if score > cos_score_cutoff:
                            # Use the original index to update df_updated
                            original_idx = group.loc[i, 'index']
                            df_updated.loc[original_idx, 'selected'] = False
                            df_updated.loc[original_idx, 'discard_reason'] = 'almost_identical_spectrum'
                            break

    return df_updated





def compare_spectra_cos(spec1, spec2, ms2_tol_da=0.02):
    """
    Compare two spectra using cosine similarity
    """
    # Create Spectrum objects
    spectrum1 = Spectrum(mz=np.array(spec1['mz_ls']),
                         intensities=np.array(spec1['intensity_ls']),
                         metadata={'precursor_mz': spec1['PEPMASS']})
    spectrum2 = Spectrum(mz=np.array(spec2['mz_ls']),
                         intensities=np.array(spec2['intensity_ls']),
                         metadata={'precursor_mz': spec2['PEPMASS']})

    # Calculate cosine similarity
    cos_greedy = CosineGreedy(tolerance=ms2_tol_da)
    score = cos_greedy.pair(spectrum1, spectrum2)

    return score['score']


def write_to_mgf(df, out_mgf):
    """
    Write the selected library to a new mgf file
    """
    df = df[df['selected']].reset_index(drop=True)
    # remove cols
    df = df.drop(['selected', 'OTHER_MATCHED_COMPOUNDS', 'OTHER_MATCHED_COMPOUNDS_NAMES', 'SPECTYPE'], axis=1)
    # move cols mz_ls and intensity_ls to the end
    df = df[[col for col in df.columns if col not in ['Num peaks', 'mz_ls', 'intensity_ls']] + ['Num peaks', 'mz_ls',
                                                                                                'intensity_ls']]

    with open(out_mgf, 'w') as file:
        # scan_idx = 1
        for idx, row in df.iterrows():
            file.write('BEGIN IONS\n')
            for key, value in row.items():
                if key == 'mz_ls':
                    for mz, intensity in zip(row['mz_ls'], row['intensity_ls']):
                        file.write(f'{mz} {intensity}\n')
                elif key == 'intensity_ls':
                    continue
                # elif key == 'SCANS':
                #     file.write(f'SCANS={scan_idx}\n')
                # elif key == 'FEATURE_ID':
                #     file.write(f'FEATURE_ID={value}\n')
                else:
                    file.write(f'{key}={value}\n')
            file.write('END IONS\n\n')
            # scan_idx += 1


def write_tsv(df, library_tsv, out_tsv):
    """
    Write the selected library to a new tsv file, for library generation on GNPS
    """
    file_path = os.path.join('input', library_tsv)
    lib_tsv = pd.read_csv(file_path, sep='\t')

    # all SCANs selected
    selected_scans = df[df['selected']]['SCANS'].tolist()

    # reserve lib_tsv with selected scans in 'EXTRACTSCAN' column
    lib_tsv = lib_tsv[lib_tsv['EXTRACTSCAN'].astype(str).isin(selected_scans)]

    # # fill 'EXTRACTSCAN' column with 1 to length of lib_tsv
    # lib_tsv['EXTRACTSCAN'] = range(1, lib_tsv.shape[0] + 1)

    # FILENAME
    lib_tsv['FILENAME'] = lib_tsv['FILENAME'].apply(lambda x: x.split('.mgf')[0] + '_filtered.mgf')

    lib_tsv.to_csv(out_tsv, sep='\t', index=False, na_rep='N/A')


def precursor_check(df, cmpd_df_dict, ms2_tol_da, prec_intensity_cutoff=10):
    """
    precursor existence check
    this is to remove the spectra that are labeled incorrectly (one spectrum being selected multiple times)
    """

    for idx, row in df.iterrows():
        if not row['selected']:
            continue
        if not row['conjugate']:
            continue
        if row['OTHER_MATCHED_COMPOUNDS'] is None:
            continue

        # get the precursor mz M+H of the two compounds
        cmpd_1_prec_mz = cmpd_df_dict.get(row['NAME_1'])
        cmpd_2_prec_mz = cmpd_df_dict.get(row['NAME_2'])
        cmpd_1_prec_int = 0
        cmpd_2_prec_int = 0
        if cmpd_1_prec_mz:
            mz_idx = np.argmin(np.abs(np.array(row['mz_ls']) - cmpd_1_prec_mz))
            if np.abs(row['mz_ls'][mz_idx] - cmpd_1_prec_mz) <= ms2_tol_da:
                cmpd_1_prec_int = row['intensity_ls'][mz_idx]
                df.loc[idx, 'cmpd_1_prec_int'] = cmpd_1_prec_int
        if cmpd_2_prec_mz:
            mz_idx = np.argmin(np.abs(np.array(row['mz_ls']) - cmpd_2_prec_mz))
            if np.abs(row['mz_ls'][mz_idx] - cmpd_2_prec_mz) <= ms2_tol_da:
                cmpd_2_prec_int = row['intensity_ls'][mz_idx]
                df.loc[idx, 'cmpd_2_prec_int'] = cmpd_2_prec_int

        df.loc[idx, 'max_prec_int'] = max(cmpd_1_prec_int, cmpd_2_prec_int)

    # sort df by _RT and _PEPMASS
    df = df.sort_values(['_RT', '_PEPMASS'])
    df['_PEPMASS'] = df['_PEPMASS'].astype(str)

    _df = df[df['conjugate']]
    discarded_scan_ls = []
    for binned_rt in _df['_RT'].unique():
        subdf = _df[_df['_RT'] == binned_rt]
        # choose a subdf with the same precursor mz
        for prec_mz in subdf['_PEPMASS'].unique():
            _subdf = subdf[subdf['_PEPMASS'] == prec_mz]

            if len(_subdf) == 1:
                continue

            # different compounds
            if len(_subdf['_NAME'].unique()) == 1:
                continue

            # maximum precursor intensity
            if _subdf['max_prec_int'].max() < prec_intensity_cutoff:
                continue

            # if multiple, sort by max_prec_int
            _subdf = _subdf.sort_values('max_prec_int', ascending=False)

            to_discard = _subdf[_subdf['max_prec_int'] == 0]['db_idx'].tolist()
            discarded_scan_ls.extend(to_discard)

    return discarded_scan_ls


def core_adduct_check(df, core_adduct_ls=None,
                      ms2_tol_da=0.05, modcos_score_cutoff=0.6, modcos_peak_cutoff=4):
    df = df[df['selected']].copy()
    # get the unique molecules
    unique_smiles = df['SMILES'].unique().tolist()
    # list of selected scan numbers
    discarded_scan_ls = []
    cmpd_adduct_to_be_removed_dict = {}  # key: RT, value: cmpd_adduct_to_be_removed
    # filter the library by (SMILES, _RT)
    for smiles in unique_smiles:
        sub_df = df[(df['SMILES'] == smiles)]
        for binned_rt in sub_df['_RT'].unique():
            _subdf = sub_df[sub_df['_RT'] == binned_rt]

            # check the subdf
            scan_ls, cmpd_adduct_to_be_removed = subdf_check(_subdf, core_adduct_ls=core_adduct_ls,
                                                             ms2_tol_da=ms2_tol_da,
                                                             modcos_score_cutoff=modcos_score_cutoff,
                                                             modcos_peak_cutoff=modcos_peak_cutoff)
            discarded_scan_ls.extend(scan_ls)
            if cmpd_adduct_to_be_removed:
                cmpd_adduct_to_be_removed_dict[binned_rt] = cmpd_adduct_to_be_removed

    return discarded_scan_ls, cmpd_adduct_to_be_removed_dict


def subdf_check(df, core_adduct_ls=None,
                ms2_tol_da=0.05, modcos_score_cutoff=0.6, modcos_peak_cutoff=4):
    """
    check the dataframe for one molecule
    :return: a list of selected scan numbers, a list of compound:adduct to be removed
    """
    if core_adduct_ls is None:
        core_adduct_ls = ['M+H', 'M-H2O+H', 'M+NH4', 'M-2H2O+H', 'M-H2O+NH4', 'M-2H2O+NH4',
                          '2M+H', '2M-H2O+H', '2M+NH4', '2M-2H2O+H', '2M-H2O+NH4', '2M-2H2O+NH4',
                          'M-H', 'M+2H', 'M-H2O+2H', 'M-2H2O+2H']

    # create a ModifiedCosine object
    modified_cosine = ModifiedCosine(tolerance=ms2_tol_da)

    subdf = df.copy()
    # check the core adducts, if any adduct in the core_adduct_ls is in the library
    subdf['core_adduct'] = subdf['_ADDUCT'].map(lambda x: x in core_adduct_ls)

    # if no core adducts are found, return empty list
    if not subdf['core_adduct'].any():
        # print('No core adducts found: ', subdf['NAME'].iloc[0], '   Found adducts: ',
        #       subdf['_ADDUCT'].unique().tolist())
        return subdf['db_idx'].tolist(), []

    # select the spectra with core adducts
    core_df = subdf[subdf['core_adduct']].copy()
    # sort by core_adduct_ls
    core_df['core_adduct'] = core_df['_ADDUCT'].map(lambda x: core_adduct_ls.index(x))
    core_df = core_df.sort_values('core_adduct')
    # get all core spectra
    core_spec_ls = []
    for idx, row in core_df.iterrows():
        core_spec = Spectrum(mz=np.array(row['mz_ls']),
                             intensities=np.array(row['intensity_ls']),
                             metadata={'precursor_mz': row['PEPMASS'],
                                       'adduct': row['_ADDUCT']})
        core_spec_ls.append(core_spec)

    unique_adducts = subdf['_ADDUCT'].unique().tolist()
    # check the adduct dependency
    if len(unique_adducts) > 1:
        allowed_adducts = unique_adducts.copy()
        if 'M-H2O+Na' in unique_adducts and 'M-H2O+H' not in unique_adducts:
            allowed_adducts.remove('M-H2O+Na')
            if 'M-2H2O+Na' in allowed_adducts:
                allowed_adducts.remove('M-2H2O+Na')
            if 'M-3H2O+Na' in allowed_adducts:
                allowed_adducts.remove('M-3H2O+Na')
        if 'M-H2O+K' in unique_adducts and 'M-H2O+H' not in unique_adducts:
            allowed_adducts.remove('M-H2O+K')
            if 'M-2H2O+K' in allowed_adducts:
                allowed_adducts.remove('M-2H2O+K')
            if 'M-3H2O+K' in allowed_adducts:
                allowed_adducts.remove('M-3H2O+K')
        if 'M-2H2O+Na' in unique_adducts and 'M-H2O+Na' not in unique_adducts and 'M-2H2O+H' in unique_adducts:
            allowed_adducts.remove('M-2H2O+Na')
            if 'M-3H2O+Na' in allowed_adducts:
                allowed_adducts.remove('M-3H2O+Na')
        if 'M-2H2O+K' in unique_adducts and 'M-H2O+K' not in unique_adducts and 'M-2H2O+H' in unique_adducts:
            allowed_adducts.remove('M-2H2O+K')
            if 'M-3H2O+K' in allowed_adducts:
                allowed_adducts.remove('M-3H2O+K')
        if 'M-3H2O+Na' in unique_adducts and 'M-2H2O+Na' not in unique_adducts and 'M-3H2O+H' in unique_adducts:
            allowed_adducts.remove('M-3H2O+Na')
        if 'M-3H2O+K' in unique_adducts and 'M-2H2O+K' not in unique_adducts and 'M-3H2O+H' in unique_adducts:
            allowed_adducts.remove('M-3H2O+K')
        if 'M-2H2O+H' in unique_adducts and 'M-H2O+H' not in unique_adducts and 'M+H' in unique_adducts:
            allowed_adducts.remove('M-2H2O+H')
            if 'M-3H2O+H' in allowed_adducts:
                allowed_adducts.remove('M-3H2O+H')
        if 'M-3H2O+H' in unique_adducts and 'M-2H2O+H' not in unique_adducts and 'M+H' in unique_adducts:
            allowed_adducts.remove('M-3H2O+H')

        if 'M-2H2O+2H' in unique_adducts and 'M-H2O+2H' not in unique_adducts and 'M+2H' in unique_adducts:
            allowed_adducts.remove('M-2H2O+2H')
            if 'M-3H2O+2H' in allowed_adducts:
                allowed_adducts.remove('M-3H2O+2H')
        if 'M-3H2O+2H' in unique_adducts and 'M-2H2O+2H' not in unique_adducts and 'M+2H' in unique_adducts:
            allowed_adducts.remove('M-3H2O+2H')

        if 'M-2H2O+NH4' in unique_adducts and 'M-H2O+NH4' not in unique_adducts and 'M-2H2O+H' not in unique_adducts:
            allowed_adducts.remove('M-2H2O+NH4')
            if 'M-3H2O+NH4' in allowed_adducts:
                allowed_adducts.remove('M-3H2O+NH4')
        if 'M-3H2O+NH4' in unique_adducts and 'M-2H2O+NH4' not in unique_adducts and 'M-3H2O+H' not in unique_adducts:
            allowed_adducts.remove('M-3H2O+NH4')

        if '2M-H2O+Na' in unique_adducts and '2M-H2O+H' not in unique_adducts:
            allowed_adducts.remove('2M-H2O+Na')
            if '2M-2H2O+Na' in allowed_adducts:
                allowed_adducts.remove('2M-2H2O+Na')
            if '2M-3H2O+Na' in allowed_adducts:
                allowed_adducts.remove('2M-3H2O+Na')
        if '2M-H2O+K' in unique_adducts and '2M-H2O+H' not in unique_adducts:
            allowed_adducts.remove('2M-H2O+K')
            if '2M-2H2O+K' in allowed_adducts:
                allowed_adducts.remove('2M-2H2O+K')
            if '2M-3H2O+K' in allowed_adducts:
                allowed_adducts.remove('2M-3H2O+K')
        if '2M-2H2O+Na' in unique_adducts and '2M-H2O+Na' not in unique_adducts and '2M-2H2O+H' in unique_adducts:
            allowed_adducts.remove('2M-2H2O+Na')
            if '2M-3H2O+Na' in allowed_adducts:
                allowed_adducts.remove('2M-3H2O+Na')
        if '2M-2H2O+K' in unique_adducts and '2M-H2O+K' not in unique_adducts and '2M-2H2O+H' in unique_adducts:
            allowed_adducts.remove('2M-2H2O+K')
            if '2M-3H2O+K' in allowed_adducts:
                allowed_adducts.remove('2M-3H2O+K')
        if '2M-3H2O+Na' in unique_adducts and '2M-2H2O+Na' not in unique_adducts and '2M-3H2O+H' in unique_adducts:
            allowed_adducts.remove('2M-3H2O+Na')
        if '2M-3H2O+K' in unique_adducts and '2M-2H2O+K' not in unique_adducts and '2M-3H2O+H' in unique_adducts:
            allowed_adducts.remove('2M-3H2O+K')
        if '2M-2H2O+H' in unique_adducts and '2M-H2O+H' not in unique_adducts and '2M+H' in unique_adducts:
            allowed_adducts.remove('2M-2H2O+H')
            if '2M-3H2O+H' in allowed_adducts:
                allowed_adducts.remove('2M-3H2O+H')
        if '2M-3H2O+H' in unique_adducts and '2M-2H2O+H' not in unique_adducts and '2M+H' in unique_adducts:
            allowed_adducts.remove('2M-3H2O+H')

        if '2M-2H2O+2H' in unique_adducts and '2M-H2O+2H' not in unique_adducts and '2M+2H' in unique_adducts:
            allowed_adducts.remove('2M-2H2O+2H')
            if '2M-3H2O+2H' in allowed_adducts:
                allowed_adducts.remove('2M-3H2O+2H')
        if '2M-3H2O+2H' in unique_adducts and '2M-2H2O+2H' not in unique_adducts and '2M+2H' in unique_adducts:
            allowed_adducts.remove('2M-3H2O+2H')

        if '2M-2H2O+NH4' in unique_adducts and '2M-H2O+NH4' not in unique_adducts and '2M-2H2O+H' not in unique_adducts:
            allowed_adducts.remove('2M-2H2O+NH4')
            if '2M-3H2O+NH4' in allowed_adducts:
                allowed_adducts.remove('2M-3H2O+NH4')
        if '2M-3H2O+NH4' in unique_adducts and '2M-2H2O+NH4' not in unique_adducts and '2M-3H2O+H' not in unique_adducts:
            allowed_adducts.remove('2M-3H2O+NH4')

        subdf['selected'] = subdf['_ADDUCT'].map(lambda x: x in allowed_adducts)

        # for unselected spectra, calculate the scores for different adducts
        for idx, row in subdf.iterrows():
            if not row['selected']:
                this_spec = Spectrum(mz=np.array(row['mz_ls']),
                                     intensities=np.array(row['intensity_ls']),
                                     metadata={'precursor_mz': row['PEPMASS'],
                                               'adduct': row['_ADDUCT']})
                # calculate the ModifiedCosine score
                for core_spec in core_spec_ls:
                    if core_spec.get('adduct') == this_spec.get('adduct'):
                        continue
                    score = modified_cosine.pair(core_spec, this_spec)
                    if score['score'] >= modcos_score_cutoff or score['matches'] >= modcos_peak_cutoff:
                        subdf.loc[idx, 'selected'] = True
                        break

        # use spectrum with unambiguous compound to filter out the ambiguous ones
        cmpd_adduct_to_be_removed = None
        _subdf = subdf[(subdf['selected']) & (subdf['OTHER_MATCHED_COMPOUNDS'].isna())]
        if len(_subdf) > 0:
            cmpd_adduct_to_be_removed = subdf['OTHER_MATCHED_COMPOUNDS_NAMES'].tolist()
            # remove nan
            cmpd_adduct_to_be_removed = [x for x in cmpd_adduct_to_be_removed if str(x) != 'nan']
            # remove empty string
            cmpd_adduct_to_be_removed = [x for x in cmpd_adduct_to_be_removed if x]
            # split by ';' and flatten the list
            cmpd_adduct_to_be_removed = [x.split(';') for x in cmpd_adduct_to_be_removed]
            cmpd_adduct_to_be_removed = [item for sublist in cmpd_adduct_to_be_removed for item in sublist]

        # db_idx of the discarded spectra
        discarded_scan_no_ls = subdf[~subdf['selected']]['db_idx'].tolist()
        return discarded_scan_no_ls, cmpd_adduct_to_be_removed
    else:
        return [], None


def isobaric_mass_adduct_check(df, cmpd_adduct_to_be_removed_dict):
    """
    isobaric mass filter, based on unique adduct existence
    """
    if cmpd_adduct_to_be_removed_dict:
        for binned_rt, cmpd_adduct_to_be_removed in cmpd_adduct_to_be_removed_dict.items():
            for cmpd_adduct in cmpd_adduct_to_be_removed:
                name, adduct, _ = cmpd_adduct.split(':')
                row_idx = df[
                    (df['_RT'] == binned_rt) & (df['NAME'] == name.strip()) & (df['ADDUCT'] == adduct.strip())].index
                if len(row_idx) > 0:
                    df.loc[row_idx, 'selected'] = False
                    df.loc[row_idx, 'discard_reason'] = 'isobaric_mass_adduct_check'
    return df


def preprocess_cmpd_df(library_csv):
    """
    Preprocess the compound dataframe, return a dictionary of compound name to precursor mass (M+H)
    """
    file_path = os.path.join('input', library_csv)
    cmpd_df = pd.read_csv(file_path)

    # rename column 'SMILES' to 'smiles'
    if 'SMILES' in cmpd_df.columns:
        cmpd_df = cmpd_df.rename(columns={'SMILES': 'smiles'})

    # check if 'formula' column exists
    if 'formula' not in cmpd_df.columns:
        print('formula column missing, converting SMILES to formula')
        cmpd_df['formula'] = cmpd_df['smiles'].apply(smiles_to_formula)

    cmpd_df['prec_mz'] = cmpd_df['formula'].apply(calc_prec_mz)

    # dictionary of mapping name to precursor mass
    cmpd_df_dict = cmpd_df.set_index('compound_name')['prec_mz'].to_dict()

    return cmpd_df_dict


def smiles_to_formula(smiles):
    # Convert SMILES to a molecule object
    mol = Chem.MolFromSmiles(smiles)

    if mol:
        # Get the molecular formula
        formula = rdMolDescriptors.CalcMolFormula(mol)
        return formula
    else:
        return None


def calc_prec_mz(formula):
    """
    Calculate the precursor mass (M+H) for a given formula string
    """
    try:
        f = Formula(formula)
        return f.monoisotopic_mass + 1.007276
    except:
        return None


def main(name_sep='_'):
    # create the output folder
    os.makedirs('output', exist_ok=True)

    metadata_folder = os.path.join('output', 'library_metadata')
    os.makedirs(metadata_folder, exist_ok=True)

    # list all the files in the input folder
    all_files = os.listdir('input')
    # filter the mgf files
    library_mgfs = [x for x in all_files if x.endswith('.mgf')]

    for library_mgf in library_mgfs:
        print(f'Processing {library_mgf}')

        base_name = library_mgf.split('.mgf')[0]

        # check the existence of the .tsv
        library_tsv = base_name + '.tsv'
        file_path = os.path.join('input', library_tsv)
        if not os.path.exists(file_path):
            print(f'tsv file missing: {library_tsv} does not exist')
            continue

        # check the existence of the .csv
        library_csv = base_name + '.csv'
        file_path = os.path.join('input', library_csv)
        if not os.path.exists(file_path):
            print(f'csv file missing: {library_csv} does not exist')
            continue
        cmpd_df_dict = preprocess_cmpd_df(library_csv)

        # main process
        df = generate_library_df(library_mgf, name_sep)
        df = select_library(df, cmpd_df_dict, base_name,
                            rt_tol=2, ms2_tol_da=0.02, prec_intensity_cutoff=10,
                            modcos_score_cutoff=0.6, modcos_peak_cutoff=4,
                            cos_score_cutoff=0.95,
                            write_df=True)

        out_mgf = os.path.join('output', base_name + '_filtered.mgf')
        write_to_mgf(df, out_mgf)
        out_tsv = os.path.join('output', base_name + '_filtered.tsv')
        write_tsv(df, library_tsv, out_tsv)


if __name__ == '__main__':
    main()
