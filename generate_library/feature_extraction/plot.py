from .raw_data_utils import MSData
from .main import init_config


def plot_all_ms2(df, file_path=None, out_dir=None):
    """
    Plot MS2 spectrum
    """
    config = init_config(mass_detect_int_tol=0)
    d = MSData()
    d.read_raw_data(file_path, config)

