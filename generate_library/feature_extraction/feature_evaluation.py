import numpy as np


def calculate_noise_level(y, intensity_threshold=0.1):
    """
    Calculate the noise level of the peak shape

    Parameters
    ----------
    y: numpy array
        Intensity
    
    Returns
    -------
    float
        noise level
    """

    y = y[y > np.max(y) * intensity_threshold]

    if len(y) < 5:
        return 0.0

    diff = np.diff(y)
    signs = np.sign(diff)
    counter = -1
    for i in range(1, len(diff)):
        if signs[i] != signs[i-1]:
            counter += 1
    if counter == -1:
        return 0.0
    
    return counter / (len(y)-2)
