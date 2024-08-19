
import json
from importlib.metadata import version


class Params:

    def __init__(self):
        """
        Function to initiate Params.
        """

        # The project
        self.project_dir = None  # Project directory, character string
        self.sample_names = None  # Absolute paths to the raw files, without extension, list of character strings
        self.sample_groups = None  # Sample groups, list of character strings
        self.sample_group_num = None  # Number of sample groups, integer
        self.sample_dir = None  # Directory for the sample information, character string
        self.single_file_dir = None  # Directory for the single file output, character string

        # MS data acquisition
        self.rt_range = [0.0, 1000.0]  # RT range in minutes, list of two floats
        self.ion_mode = "positive"  # Ionization mode, "positive" or "negative", character string

        # Feature detection
        self.mz_tol_ms1 = 0.01  # m/z tolerance for MS1, default is 0.01
        self.mz_tol_ms2 = 0.015  # m/z tolerance for MS2, default is 0.015
        self.int_tol = 30000  # Intensity tolerance, default is 30000 for Orbitrap and 1000 for other instruments, integer
        self.roi_gap = 30  # Gap within a feature, default is 30 (i.e. 30 consecutive scans without signal), integer
        self.ppr = 0.8  # Peak peak correlation threshold for feature grouping, default is 0.8

    def set_default(self, ms_type, ion_mode):
        """
        Set the parameters by the type of MS.
        --------------------------------------
        ms_type: character string
            The type of MS, "orbitrap" or "tof".
        ion_mode: character string
            The ionization mode, "positive" or "negative".
        """

        if ms_type == "orbitrap":
            self.int_tol = 30000
        elif ms_type == "tof":
            self.int_tol = 1000
        if ion_mode == "positive":
            self.ion_mode = "positive"
        elif ion_mode == "negative":
            self.ion_mode = "negative"

    def output_parameters(self, path, format="json"):
        """
        Output the parameters to a file.
        ---------------------------------

        Parameters
        ----------
        path : str
            The path to the output file.
        format : str
            The format of the output file. "json" is only supported for now.
        """

        if format == "json":
            parameters = {}
            # obtain the version of the package
            parameters["MassCube_version"] = version("masscube")

            for key, value in self.__dict__.items():
                if key != "project_dir":
                    parameters[key] = value
            with open(path, 'w') as f:
                json.dump(parameters, f)
        else:
            raise ValueError("The output format is not supported.")


def find_ms_info(file_name):
    """
    Find the type of MS, ion mode, centroided from the raw file.
    """

    ms_type = 'tof'
    ion_mode = 'positive'
    centroid = False

    with open(file_name, 'r') as f:
        for i, line in enumerate(f):
            if 'orbitrap' in line.lower():
                ms_type = 'orbitrap'
            if 'negative' in line.lower():
                ion_mode = 'negative'
            if "centroid spectrum" in line.lower() or 'centroided="1"' in line.lower():
                centroid = True
            if i > 200:
                break

    return ms_type, ion_mode, centroid
