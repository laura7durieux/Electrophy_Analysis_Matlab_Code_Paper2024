#### main script #####

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Source.DLCResult import DLCResult, Column, Head
from Source.mvt_sorting import distance_to_movingness, average_movingness_to_movingness, coord_to_distance
from Source.chart import display_chart, display_double_chart, display_threshold_chart

# Loading the data


def analyze_mvt_sleep(path, name_saved):
    dlc = DLCResult(path)
    # manual inputs
    num_frame_per_sec = 30
    path_results = 'results_sleep\\'

    # Analyse dist for Hat as example
    data = dlc.extract_coords(Head.Hat)
    dist_hat = coord_to_distance(data)

    print(len(dist_hat))
    # calculate time
    total_time = np.size(dist_hat, 0) / num_frame_per_sec
    time = np.arange(0, total_time, 1 / num_frame_per_sec)  # in sec
    time = time.reshape((-1, 1)) # transpose

    # init variable for calculation of all bodyparts
    row_count = np.size(dist_hat, 0)
    column_count = len(Head)

    dist = np.zeros((np.size(dist_hat,0), 2)) # based on the time length

    movingness_per_column = np.zeros((row_count, column_count))
    likelihood_per_column = np.zeros((row_count, column_count))

    # Process for all point (bodyparts)
    for point in Head:
        data = dlc.extract_coords(point)
        dist = coord_to_distance(data)

        # Calculate the threshold (based on the 98.5 percentile)
        quartile80 = np.percentile(dist[:,0], 98.5)

        threshold = quartile80
        print('threshold=', threshold)

        # Calculate the binary code 0 don't move 1 moving based on the threshold
        movingness = distance_to_movingness(dist, threshold)

        # # Display plot with the threshold
        # display_threshold_chart(dist, time, threshold, path_results +"sleep_Treshold_" + str(point))
        #
        # # save the  chart representing the pixel distance and the binary code for each bodyparts
        # display_double_chart(dist[:,0], movingness[:, 0],time, path_results+"Binary_sleep_"+str(point))

        # Stock the result into an array (time x type)
        movingness_per_column[:, int(point)] = movingness[:, 0]
        likelihood_per_column[:, int(point)] = movingness[:, 1]

    # Calculate the average of the bodyparts weighted by the likelihood of each points
    average_movingness = np.zeros((row_count, 1))

    for i in range(0, row_count):
        # print('movingness_per_column=', movingness_per_column[i])
        # print('likelihood_per_column=', likelihood_per_column[i])
        if np.sum(likelihood_per_column[i]) > 0:
            average_movingness[i, 0] = np.average(movingness_per_column[i], weights=likelihood_per_column[i])
        else:
            average_movingness[i, 0] = 0.5  # Not sur if it's moving or not !

    print(average_movingness)

    # do the threshold at 0.5 on smoothed data
    final_movingness = average_movingness_to_movingness(average_movingness, 0.5)

    # Do the plots
    # display_double_chart(average_movingness, final_movingness, time, path_results + "Sleep _ Average on the body part after and before the threshold")
    # display_double_chart(dist_hat, final_movingness,time, path_results + "Sleep _ Threshold average compare to Hat")
    if len(final_movingness) != len(time):
        time = time[0:-1]

    # transform array to data frame
    movingness_sleep_output = pd.DataFrame(np.transpose(np.concatenate((final_movingness, time), axis=1)))

    print('movingness_sleep_output'+str(movingness_sleep_output))

    movingness_sleep_output.to_csv(path_results + name_saved + '_sleep-movingness.csv', index = False, header=False)


# Do a list of all vids for results
import os
path_vids = "D:\\Documents\\My_GitHub\\Behavior_analysis\\Data\\Elect4\\"

list_vids_full = []
name_saved = []
i = 1
for f in os.listdir(path_vids):
    name, ext = os.path.splitext(f)
    if ext == '.csv':
        list_vids_full.append(path_vids+f)
        name_saved.append(name)
        i = i+1


count = 0
for vids in list_vids_full:
    print(vids)
    print(name_saved[count])
    analyze_mvt_sleep(vids, name_saved[count])
    count = count + 1
