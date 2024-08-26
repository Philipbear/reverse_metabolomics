from pyteomics import mzml, mzxml
import numpy as np
import os
import pandas as pd

from .config import Params
from .peak_detection import find_rois, cut_roi


class MSData:
    """
    A class that models a single file (mzML or mzXML) and processes the raw data.
    """

    def __init__(self):
        """
        Function to initiate MSData.
        ----------------------------------------------------------
        """

        self.scans = []  # A list of Scan objects
        self.ms1_rt_seq = []  # Retention times of all MS1 scans
        self.bpc_int = []  # Intensity of the BPC
        self.rois = []  # A list of ROIs
        self.params = None  # A Params object
        self.file_name = None  # File name of the raw data without extension
        self.roi_mz_seq = None  # m/z of all ROIs
        self.roi_rt_seq = None  # Retention time of all ROIs

    def read_raw_data(self, file_name, params, read_ms2=True, centroid=True):
        """
        Function to read raw data to MS1 and MS2 (if available)
        (supported by pyteomics package).

        Parameters
        ----------------------------------------------------------
        file_name: str
            File name of raw MS data (mzML or mzXML).
        params: Params object
            A Params object that contains the parameters.
        """

        self.params = params

        if os.path.isfile(file_name):
            # get extension from file name
            ext = os.path.splitext(file_name)[1]

            if ext.lower() != ".mzml" and ext.lower() != ".mzxml":
                raise ValueError("Unsupported raw data format. Raw data must be in mzML or mzXML.")

            self.file_name = os.path.splitext(os.path.basename(file_name))[0]

            if ext.lower() == ".mzml":
                with mzml.MzML(file_name) as reader:
                    self.extract_scan_mzml(reader, int_tol=params.int_tol, read_ms2=read_ms2, centroid=centroid)
            elif ext.lower() == ".mzxml":
                with mzxml.MzXML(file_name) as reader:
                    self.extract_scan_mzxml(reader, int_tol=params.int_tol, read_ms2=read_ms2, centroid=centroid)
        else:
            print("File does not exist.")

    def extract_scan_mzml(self, spectra, int_tol=0, read_ms2=True, centroid=True):
        """
        Function to extract all scans and convert them to Scan objects.

        Parameters
        ----------------------------------------------------------
        spectra: pyteomics object
            An iteratable object that contains all MS1 and MS2 scans.
        """

        idx = 0  # Scan number
        self.ms1_idx = []  # MS1 scan index
        self.ms2_idx = []  # MS2 scan index

        rt_unit = spectra[0]['scanList']['scan'][0]['scan start time'].unit_info

        # Iterate over all scans
        for spec in spectra:
            # Get the retention time and convert to minute
            try:
                rt = spec['scanList']['scan'][0]['scan start time']
            except:
                rt = spec['scanList']['scan'][0]['scan time']

            if rt_unit == 'second':
                rt = rt / 60

            # Check if the retention time is within the range
            if self.params.rt_range[0] < rt < self.params.rt_range[1]:
                if spec['ms level'] == 1:
                    temp_scan = Scan(level=1, scan=idx, rt=rt)
                    mz_array = np.array(spec['m/z array'], dtype=np.float64)
                    int_array = np.array(spec['intensity array'], dtype=np.int64)
                    mz_array = mz_array[int_array > int_tol]
                    int_array = int_array[int_array > int_tol]

                    if len(mz_array) == 0:
                        continue

                    if centroid:
                        mz_array, int_array = _centroid(mz_array, int_array)

                    temp_scan.add_info_by_level(mz_seq=mz_array, int_seq=int_array)
                    self.ms1_idx.append(idx)

                    # update base peak chromatogram
                    self.bpc_int.append(np.max(int_array))
                    self.ms1_rt_seq.append(rt)

                elif spec['ms level'] == 2 and read_ms2:
                    temp_scan = Scan(level=2, scan=idx, rt=rt)
                    precursor_mz = spec['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0][
                        'selected ion m/z']
                    peaks = np.array([spec['m/z array'], spec['intensity array']], dtype=np.float64).T
                    temp_scan.add_info_by_level(precursor_mz=precursor_mz, peaks=peaks)

                    self.ms2_idx.append(idx)

                self.scans.append(temp_scan)
                idx += 1

        self.ms1_rt_seq = np.array(self.ms1_rt_seq)

    def extract_scan_mzxml(self, spectra, int_tol, read_ms2=True, centroid=True):
        """
        Function to extract all scans and convert them to Scan objects.

        Parameters
        ----------------------------------------------------------
        spectra: pyteomics object
            An iteratable object that contains all MS1 and MS2 scans.
        """

        idx = 0  # Scan number
        self.ms1_idx = []  # MS1 scan index
        self.ms2_idx = []  # MS2 scan index

        rt_unit = spectra[0]["retentionTime"].unit_info

        # Iterate over all scans
        for spec in spectra:
            # Get the retention time and convert to minute
            rt = spec["retentionTime"]  # retention time of mzXML is in minute

            if rt_unit == 'second':
                rt = rt / 60

            # Check if the retention time is within the range
            if self.params.rt_range[0] < rt < self.params.rt_range[1]:
                if spec['msLevel'] == 1:
                    temp_scan = Scan(level=1, scan=idx, rt=rt)
                    mz_array = np.array(spec['m/z array'], dtype=np.float64)
                    int_array = np.array(spec['intensity array'], dtype=np.int64)

                    mz_array = mz_array[int_array > int_tol]
                    int_array = int_array[int_array > int_tol]

                    if centroid:
                        mz_array, int_array = _centroid(mz_array, int_array)

                    temp_scan.add_info_by_level(mz_seq=mz_array, int_seq=int_array)
                    self.ms1_idx.append(idx)

                    # update base peak chromatogram
                    self.bpc_int.append(np.max(int_array))
                    self.ms1_rt_seq.append(rt)

                elif spec['msLevel'] == 2 and read_ms2:
                    temp_scan = Scan(level=2, scan=idx, rt=rt)
                    precursor_mz = spec['precursorMz'][0]['precursorMz']
                    peaks = np.array([spec['m/z array'], spec['intensity array']], dtype=np.float64).T
                    temp_scan.add_info_by_level(precursor_mz=precursor_mz, peaks=peaks)

                    self.ms2_idx.append(idx)

                self.scans.append(temp_scan)
                idx += 1

        self.ms1_rt_seq = np.array(self.ms1_rt_seq)

    def drop_ion_by_int(self):
        """
        Function to drop ions by intensity.

        Parameters
        ----------------------------------------------------------
        tol: int
            Intensity tolerance.
        """

        for idx in self.ms1_idx:
            self.scans[idx].mz_seq = self.scans[idx].mz_seq[self.scans[idx].int_seq > self.params.int_tol]
            self.scans[idx].int_seq = self.scans[idx].int_seq[self.scans[idx].int_seq > self.params.int_tol]

    def find_rois(self):
        """
        Function to find ROI in MS1 scans.

        Parameters
        ----------------------------------------------------------
        params: Params object
            A Params object that contains the parameters.
        """

        self.rois = find_rois(self)

        for roi in self.rois:
            roi.sum_roi()

        self.rois.sort(key=lambda x: x.mz)

    def cut_rois(self):
        """
        Function to cut ROI into smaller pieces.
        """

        self.rois = [cut_roi(r, int_tol=self.params.int_tol) for r in self.rois]
        tmp = []
        for roi in self.rois:
            tmp.extend(roi)

        self.rois = tmp

    def summarize_roi(self):
        """
        Function to process ROIs.

        Parameters
        ----------------------------------------------------------
        params: Params object
            A Params object that contains the parameters.
        """

        for roi in self.rois:
            roi.sum_roi()

        # sort rois by m/z
        self.rois.sort(key=lambda x: x.mz)

        # index the rois
        for idx in range(len(self.rois)):
            self.rois[idx].id = idx

        # extract mz and rt of all rois for further use (feature grouping)
        self.roi_mz_seq = np.array([roi.mz for roi in self.rois])
        self.roi_rt_seq = np.array([roi.rt for roi in self.rois])

        # allocate ms2 to rois
        self.allocate_ms2_to_rois()

        # find best ms2 for each roi
        for roi in self.rois:
            if len(roi.ms2_seq) > 0:
                roi.best_ms2 = find_best_ms2(roi.ms2_seq)

    def allocate_ms2_to_rois(self):
        for i in self.ms2_idx:
            if len(self.scans[i].peaks) == 0:
                continue
            idx = np.where(np.abs(self.roi_mz_seq - self.scans[i].precursor_mz) < self.params.mz_tol_ms1)[0]
            matched_rois = []
            for j in idx:
                if self.rois[j].rt_seq[0] < self.scans[i].rt < self.rois[j].rt_seq[-1]:
                    matched_rois.append(self.rois[j])
            if len(matched_rois) > 0:
                # assign ms2 to the roi with the highest peak height
                matched_rois[np.argmax([roi.peak_height for roi in matched_rois])].ms2_seq.append(self.scans[i])

    def drop_rois_without_ms2(self):
        """
        Function to drop ROIs without MS2.
        """

        self.rois = [roi for roi in self.rois if len(roi.ms2_seq) > 0]

    def drop_rois_by_length(self, length=5):
        """
        Function to drop ROIs by length.
        """

        self.rois = [roi for roi in self.rois if roi.length >= length]

    def _discard_isotopes(self):
        """
        Function to discard isotopes.
        """

        self.rois = [roi for roi in self.rois if not roi.is_isotope]

        for idx in range(len(self.rois)):
            self.rois[idx].id = idx

    def output_single_file(self, save=True, out_dir=None):
        """
        Function to generate a report for rois.
        """
        result = []

        for roi in self.rois:

            # iso_peaks = np.column_stack((roi.isotope_mz_seq, roi.isotope_int_seq))

            temp = [roi.id, roi.mz, roi.rt, roi.length, roi.rt_seq[0],
                    roi.rt_seq[-1], roi.peak_area, roi.peak_height, roi.charge_state,
                    # roi.is_isotope, roi.isotope_id_seq, iso_peaks,
                    # roi.adduct_type,
                    # int(roi.adduct_group_no) if roi.adduct_group_no is not None else None,
                    roi.best_ms2.scan + 1 if roi.best_ms2 is not None else None,
                    roi.best_ms2.peaks if roi.best_ms2 is not None else None,
                    [ms2.scan + 1 for ms2 in roi.ms2_seq] if roi.ms2_seq else None]

            result.append(temp)

        # convert result to a pandas dataframe
        columns = ["ID", "m/z", "RT", "length", "RT_start", "RT_end", "peak_area", "peak_height",
                   "charge", #"is_isotope", "isotope_IDs", "isotopes", "adduct", "adduct_group",
                   "best_MS2_scan_idx", "MS2", "all_MS2_scan_idx"]

        df = pd.DataFrame(result, columns=columns)

        # save
        if save:
            if not out_dir:
                out_dir = self.params.single_file_dir

            os.makedirs(out_dir, exist_ok=True)
            path = os.path.join(out_dir, self.file_name + "_feature_table.tsv")
            df.to_csv(path, sep="\t", index=False)

        return df

    def get_eic_data(self, target_mz, target_rt=None, mz_tol=0.005, rt_tol=0.3, rt_range=None):
        """
        To get the EIC data of a target m/z.

        Parameters
        ----------
        target_mz: float
            Target m/z.
        mz_tol: float
            m/z tolerance.
        target_rt: float
            Target retention time.
        rt_tol: float
            Retention time tolerance.

        Returns
        -------
        eic_rt: numpy array
            Retention time of the EIC.
        eic_int: numpy array
            Intensity of the EIC.
        eic_mz: numpy array
            m/z of the EIC.
        eic_scan_idx: numpy array
            Scan index of the EIC.
        """

        eic_rt = []
        eic_int = []
        eic_mz = []
        eic_scan_idx = []

        if target_rt is not None:
            rt_range = [target_rt - rt_tol, target_rt + rt_tol]
        elif rt_range is None:
            rt_range = [0, np.inf]

        for i in self.ms1_idx:
            if self.scans[i].rt > rt_range[0] and self.scans[i].rt < rt_range[1]:
                mz_diff = np.abs(self.scans[i].mz_seq - target_mz)
                if len(mz_diff) > 0 and np.min(mz_diff) < mz_tol:
                    eic_rt.append(self.scans[i].rt)
                    eic_int.append(self.scans[i].int_seq[np.argmin(mz_diff)])
                    eic_mz.append(self.scans[i].mz_seq[np.argmin(mz_diff)])
                    eic_scan_idx.append(i)
                else:
                    eic_rt.append(self.scans[i].rt)
                    eic_int.append(0)
                    eic_mz.append(0)
                    eic_scan_idx.append(i)

            if self.scans[i].rt > rt_range[1]:
                break

        eic_rt = np.array(eic_rt)
        eic_int = np.array(eic_int)
        eic_mz = np.array(eic_mz)
        eic_scan_idx = np.array(eic_scan_idx)

        return eic_rt, eic_int, eic_mz, eic_scan_idx

    def find_ms2_by_mzrt(self, mz_target, rt_target, mz_tol=0.01, rt_tol=0.3, return_best=False):
        """
        Function to find MS2 scan by precursor m/z and retention time.

        Parameters
        ----------------------------------------------------------
        mz_target: float
            Precursor m/z.
        rt_target: float
            Retention time.
        mz_tol: float
            m/z tolerance.
        rt_tol: float
            Retention time tolerance.
        return_best: bool
            True: only return the best MS2 scan with the highest intensity.
            False: return all MS2 scans as a list.
        """

        matched_ms2 = []

        for id in self.ms2_idx:
            rt = self.scans[id].rt

            if rt < rt_target - rt_tol:
                continue

            mz = self.scans[id].precursor_mz

            if abs(mz - mz_target) < mz_tol and abs(rt - rt_target) < rt_tol:
                matched_ms2.append(self.scans[id])

            if rt > rt_target + rt_tol:
                break

        if return_best:
            if len(matched_ms2) > 1:
                total_ints = [np.sum(ms2.peaks[:, 1]) for ms2 in matched_ms2]
                return matched_ms2[np.argmax(total_ints)]
            elif len(matched_ms2) == 1:
                return matched_ms2[0]
            else:
                return None
        else:
            return matched_ms2

    def find_roi_by_mzrt(self, mz_target, rt_target=None, mz_tol=0.01, rt_tol=0.3):
        """
        Function to find roi by precursor m/z and retention time.

        Parameters
        ----------------------------------------------------------
        mz_target: float
            Precursor m/z.
        rt_target: float
            Retention time.
        mz_tol: float
            m/z tolerance.
        rt_tol: float
            Retention time tolerance.
        """

        if rt_target is None:
            tmp = np.abs(self.roi_mz_seq - mz_target) < mz_tol
            found_roi = [self.rois[i] for i in np.where(tmp)[0]]
        else:
            tmp1 = np.abs(self.roi_mz_seq - mz_target) < mz_tol
            tmp2 = np.abs(self.roi_rt_seq - rt_target) < rt_tol
            tmp = np.logical_and(tmp1, tmp2)
            found_roi = [self.rois[i] for i in np.where(tmp)[0]]

        return found_roi

    def find_ms1_scan_by_rt(self, rt_target):
        """
        Function to find a MS1 scan by retention time.

        Parameters
        ----------------------------------------------------------
        rt_target: float
            Retention time.
        """

        idx = np.argmin(np.abs(self.ms1_rt_seq - rt_target))
        return self.scans[self.ms1_idx[idx]]

    def correct_retention_time(self, f):
        """
        Function to correct retention time.

        Parameters
        ----------------------------------------------------------
        f: interp1d object
            A function to correct retention time.
        """

        all_rts = np.array([s.rt for s in self.scans])
        all_rts = f(all_rts)
        for i in range(len(self.scans)):
            self.scans[i].rt = all_rts[i]


class Scan:
    """
    A class that represents an MS scan.
    A MS1 spectrum has properties including:
        scan number, retention time,
        m/z and intensities.
    A MS2 spectrum has properties including:
        scan number, retention time,
        precursor m/z, product m/z and intensities.
    """

    def __init__(self, level=None, scan=None, rt=None):
        """
        Parameters
        ----------------------------------------------------------
        level: int
            Level of MS scan.
        scan: int
            Scan number.
        rt: float
            Retention time.
        """

        self.level = level
        self.scan = scan
        self.rt = rt

        # for MS1 scans:
        self.mz_seq = None
        self.int_seq = None

        # for MS2 scans:
        self.precursor_mz = None
        self.peaks = None

    def add_info_by_level(self, **kwargs):
        """
        Function to add scan information by level.
        """

        if self.level == 1:
            self.mz_seq = kwargs['mz_seq']
            self.int_seq = np.int64(kwargs['int_seq'])

        elif self.level == 2:
            self.precursor_mz = kwargs['precursor_mz']
            self.peaks = kwargs['peaks']

    def show_scan_info(self):
        """
        Function to print a scan's information.

        Parameters
        ----------------------------------------------------------
        scan: MS1Scan or MS2Scan object
            A MS1Scan or MS2Scan object.
        """

        print("Scan number: " + str(self.scan))
        print("Retention time: " + str(self.rt))

        if self.level == 1:
            print("m/z: " + str(np.around(self.mz_seq, decimals=4)))
            print("Intensity: " + str(np.around(self.int_seq, decimals=0)))

        elif self.level == 2:
            # keep 4 decimal places for m/z and 0 decimal place for intensity
            print("Precursor m/z: " + str(np.round(self.precursor_mz, decimals=4)))
            print(self.peaks)


def _centroid(mz_seq, int_seq, mz_tol=0.005):
    """
    Function to centroid the m/z and intensity sequences.

    Parameters
    ----------
    mz_seq: numpy array or list
        m/z sequence.
    int_seq: numpy array or list
        Intensity sequence.
    mz_tol: float
        m/z tolerance for centroiding. Default is 0.005.
    """

    diff = np.diff(mz_seq)
    tmp = np.where(diff < mz_tol)[0]

    if len(tmp) == 0:
        return mz_seq, int_seq

    mz_seq = list(mz_seq)
    int_seq = list(int_seq)
    for i in tmp[::-1]:
        mz_seq[i] = (mz_seq[i] * int_seq[i] + mz_seq[i + 1] * int_seq[i + 1]) / (int_seq[i] + int_seq[i + 1])
        int_seq[i] += int_seq[i + 1]
        mz_seq.pop(i + 1)
        int_seq.pop(i + 1)

    return np.array(mz_seq), int_seq


def read_raw_file_to_obj(file_name, params=None, int_tol=1000, centroid=True, read_ms2=True, print_summary=False):
    """
    Read a raw file to a MSData object.
    It's a useful function for data visualization or brief data analysis.

    Parameters
    ----------
    file_name : str
        The file name of the raw file.
    params : Params object
        A Params object that contains the parameters.
    int_tol : float
        Intensity tolerance for dropping ions.
    centroid : bool
        True: centroid the raw data.
        False: do not centroid the raw data.
    print_summary : bool
        True: print the summary of the raw data.
        False: do not print the summary of the raw data.

    Returns
    -------
    d : MSData object
        A MSData object.
    """

    # create a MSData object
    d = MSData()

    # read raw data
    if params is None:
        params = Params()
        params.int_tol = int_tol

    d.read_raw_data(file_name, params=params, read_ms2=read_ms2, centroid=centroid)

    if print_summary:
        print("Number of MS1 scans: " + str(len(d.ms1_idx)), "Number of MS2 scans: " + str(len(d.ms2_idx)))

    return d


def find_best_ms2(ms2_seq):
    """
    Function to find the best MS2 spectrum for a list of MS2 spectra.
    """

    if len(ms2_seq) > 0:
        total_ints = [np.sum(ms2.peaks[:, 1]) for ms2 in ms2_seq]
        if np.max(total_ints) == 0:
            return None
        else:
            return ms2_seq[max(range(len(total_ints)), key=total_ints.__getitem__)]
    else:
        return None


def write_peaks(mz_arr, intensity_arr):
    peak_str = ""
    for i in range(len(mz_arr)):
        peak_str += str(round(mz_arr[i], 4)) + "," + str(intensity_arr[i]) + ";"
    return peak_str[:-1]
