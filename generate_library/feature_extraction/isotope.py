import numpy as np


def annotate_isotope(d, mz_tol=0.015, rt_tol=0.1, valid_intensity_ratio_range=[0.001, 1.2]):
    """
    Function to annotate isotopes.
    """

    # rank the rois (d.rois) in each file by m/z
    d.rois.sort(key=lambda x: x.mz)

    for r in d.rois:
        
        if r.is_isotope:
            continue

        r.isotope_int_seq = [r.peak_height]
        r.isotope_mz_seq = [r.mz]

        last_mz = r.mz

        # check if the currect ion is double charged
        r.charge_state = 1
        target_mz = r.mz + 1.003355/2
        v = np.where(np.logical_and(np.abs(d.roi_mz_seq - target_mz) < mz_tol, np.abs(d.roi_rt_seq - r.rt) < rt_tol))[0]
        if len(v) > 0:
            # an isotope can't have intensity 1.2 fold or higher than M0 or 1% lower than the M0
            v = [v[i] for i in range(len(v)) if d.rois[v[i]].peak_height < valid_intensity_ratio_range[1] * r.peak_height]
            v = [v[i] for i in range(len(v)) if d.rois[v[i]].peak_height > 0.01*r.peak_height]
            if len(v) > 0:
                r.charge_state = 2

        isotope_id_seq = []
        target_mz = r.mz + 1.003355/r.charge_state
        total_int = r.peak_height

        # find roi using isotope list
        i = 0
        while i < 5:   # maximum 5 isotopes
            # if isotope is not found in 2.2 m/z, stop searching
            if target_mz - last_mz > 2.2:
                break

            v = np.where(np.logical_and(np.abs(d.roi_mz_seq - target_mz) < mz_tol, np.abs(d.roi_rt_seq - r.rt) < rt_tol))[0]

            if len(v) == 0:
                i += 1
                continue

            # an isotope can't have intensity 1.2 fold or higher than M0 or 0.1% lower than the M0
            v = [v[i] for i in range(len(v)) if d.rois[v[i]].peak_height < valid_intensity_ratio_range[1] * r.peak_height]
            v = [v[i] for i in range(len(v)) if d.rois[v[i]].peak_height > valid_intensity_ratio_range[0] * total_int]

            if len(v) == 0:
                i += 1
                continue
            
            # for high-resolution data, C and N isotopes can be separated and need to be summed
            total_int = np.sum([d.rois[vi].peak_height for vi in v])
            last_mz = np.mean([d.rois[vi].mz for vi in v])

            r.isotope_mz_seq.append(last_mz)
            r.isotope_int_seq.append(total_int)
            
            for vi in v:
                d.rois[vi].is_isotope = True
                isotope_id_seq.append(d.rois[vi].id)

            target_mz = last_mz + 1.003355/r.charge_state
            i += 1

        r.charge_state = get_charge_state(r.isotope_mz_seq)
        r.isotope_id_seq = isotope_id_seq


def get_charge_state(mz_seq):
    
    if len(mz_seq) < 2:
        return 1
    else:
        mass_diff = mz_seq[1] - mz_seq[0]

        # check mass diff is closer to 1 or 0.5 | note, mass_diff can be larger than 1
        if abs(mass_diff - 1) < abs(mass_diff - 0.5):
            return 1
        else:
            return 2


