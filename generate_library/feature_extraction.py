from feature_extraction import feature_extraction_single

feature_extraction_single(file_path='/Users/shipei/Documents/projects/reverse_metabolomics/reverse_metabolomics/filter_library/input/AP_68.mzML',
                          mass_detect_int_tol=5e4,
                          peak_cor_rt_tol=0.025, min_ppc=0.8,
                          ms1_tol=0.01, ms2_tol=0.015,
                          out_dir=None)
