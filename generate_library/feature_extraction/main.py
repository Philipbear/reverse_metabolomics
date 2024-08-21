import os

from .isotope import annotate_isotope
from .adduct import annotate_adduct
from .config import Params, find_ms_info
from .raw_data_utils import MSData


def feature_extraction_single(file_path, mass_detect_int_tol=5e4,
                              peak_cor_rt_tol=0.025, min_ppc=0.8,
                              ms1_tol=0.005, ms2_tol=0.015,
                              save=True, out_dir=None):
    """
    feature detection from a single .mzML file
    """

    ms_type, ion_mode, centroid = find_ms_info(file_path)

    # init a new config object
    config = init_config(ion_mode=ion_mode, mz_tol_ms1=ms1_tol, mz_tol_ms2=ms2_tol,
                         peak_cor_rt_tol=peak_cor_rt_tol, min_ppc=min_ppc,
                         mass_detect_int_tol=mass_detect_int_tol)

    # create a MSData object
    d = MSData()

    # read raw data
    print('Reading raw data...')
    d.read_raw_data(file_path, config)

    # detect region of interests (ROIs)
    print('Detecting ROIs...')
    d.find_rois()

    # cut ROIs
    d.cut_rois()

    # label short ROIs, find the best MS2, and sort ROIs by m/z
    d.summarize_roi()

    # annotate isotopes, adducts
    # print('Annotating isotopes and adducts...')
    # annotate_isotope(d)
    # annotate_adduct(d)

    # output single file to a tsv file, in the same directory as the raw file
    print('Generating feature table...')
    if out_dir is None:
        out_dir = os.path.dirname(file_path)

    df = d.output_single_file(save, out_dir)

    return df


def init_config(ion_mode="positive",
                mz_tol_ms1=0.005, mz_tol_ms2=0.015,
                peak_cor_rt_tol=0.025, min_ppc=0.8,
                mass_detect_int_tol=5e4):
    # init
    config = Params()

    ##########################
    # MS data acquisition
    config.ion_mode = ion_mode  # Ionization mode, "positive" or "negative", default is "positive"
    config.mz_tol_ms1 = mz_tol_ms1  # m/z tolerance for MS1, default is 0.01
    config.mz_tol_ms2 = mz_tol_ms2  # m/z tolerance for MS2, default is 0.015
    config.int_tol = mass_detect_int_tol  # Mass detection intensity tolerance, default is 5e4
    config.ppr = min_ppc  # Peak-peak correlation threshold for feature grouping, default is 0.8
    config.peak_cor_rt_tol = peak_cor_rt_tol  # RT tolerance for peak correlation, default is 0.025

    return config
