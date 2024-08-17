import os

from .isotope import annotate_isotope
from .adduct import annotate_adduct
from .config import Params, find_ms_info
from .raw_data_utils import MSData


def feature_extraction_single(file_path, mass_detect_int_tol=5e4,
                              peak_cor_rt_tol=0.025, min_ppc=0.8,
                              ms1_tol=0.01, ms2_tol=0.015,
                              out_dir=None):
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
    print('Annotating isotopes and adducts...')
    annotate_isotope(d)
    annotate_adduct(d)

    # output single file to a tsv file, in the same directory as the raw file
    if out_dir is None:
        out_dir = os.path.dirname(file_path)

    d.output_single_file(out_dir)

    return d


def init_config(ion_mode="positive",
                mz_tol_ms1=0.01, mz_tol_ms2=0.015,
                peak_cor_rt_tol=0.025, min_ppc=0.8,
                mass_detect_int_tol=5e4):
    # init
    config = Params()

    ##########################
    # The project
    config.project_dir = None  # Project directory, character string
    config.sample_names = None  # Absolute paths to the raw files, without extension, list of character strings
    config.sample_groups = None  # Sample groups, list of character strings
    config.sample_group_num = None  # Number of sample groups, integer
    config.sample_dir = None  # Directory for the sample information, character string
    config.single_file_dir = None  # Directory for the single file output, character string

    # MS data acquisition
    config.rt_range = [0.0, 1000.0]  # RT range in minutes, list of two floats
    config.ion_mode = ion_mode  # Ionization mode, "positive" or "negative", character string

    # Feature detection
    config.mz_tol_ms1 = mz_tol_ms1  # m/z tolerance for MS1, default is 0.01
    config.mz_tol_ms2 = mz_tol_ms2  # m/z tolerance for MS2, default is 0.015
    config.int_tol = mass_detect_int_tol  # Intensity tolerance, default is 30000 for Orbitrap and 1000 for other instruments, integer
    config.roi_gap = 30  # Gap within a feature, default is 30 (i.e. 30 consecutive scans without signal), integer
    config.ppr = min_ppc  # Peak-peak correlation threshold for feature grouping, default is 0.8
    config.peak_cor_rt_tol = peak_cor_rt_tol  # RT tolerance for peak correlation, default is 0.2

    # Parameters for output
    config.output_single_file = False  # Whether to output the processed individual files to a csv file

    return config
