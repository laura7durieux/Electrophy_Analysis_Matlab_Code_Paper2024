import numpy as np


def smoothing_by_time(to_smooth,division_factor):
    raw_to_smooth=np.size(to_smooth,0)
    length_final_to_smooth = raw_to_smooth//division_factor # floor division
    smoothed = np.zeros((length_final_to_smooth,1))


    for i in range(0,length_final_to_smooth):
        smoothed[i] = np.mean(to_smooth[i:i+division_factor])

    return smoothed