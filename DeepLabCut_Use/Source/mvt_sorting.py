import math

import numpy as np


def coord_to_distance(dataset):
    col_dataset = np.size(dataset,1)
    print(col_dataset)
    dataArray = dataset

    raw_dataset = np.size(dataArray, 0)

    # distance = racine((x2-x1)²+(y2-y1)²)
    dist = np.zeros((raw_dataset - 1, 2))
    for i in range(1, raw_dataset - 1):
        x1 = float(dataArray[i, 0])
        x2 = float(dataArray[i + 1, 0])
        y1 = float(dataArray[i, 1])
        y2 = float(dataArray[i + 1, 1])

        # On veut stocker le pourcentage de certitude avec la distance.
        # On multiplie l'un avec l'autre pour avoir quelque chose qui colle
        likelihood1 = float(dataArray[i, 2])
        likelihood2 = float(dataArray[i + 1, 2])

        distance = math.sqrt(((x2 - x1) ** 2) + ((y2 - y1) ** 2))
        dist[i - 1, 0] = distance
        dist[i - 1, 1] = likelihood1 * likelihood2
    return dist


def distance_to_movingness(dist, threshold):
    dist_size = np.size(dist, 0)

    # movingness = value > threshold
    movingness = np.zeros((dist_size, 2))
    for i in range(0, dist_size):
        movingness_for_i = 0 if dist[i, 0] < threshold else 1
        likelihood_for_i = dist[i, 1]
        movingness[i, 0] = movingness_for_i
        movingness[i, 1] = likelihood_for_i
    return movingness


def average_movingness_to_movingness(average_movingness, threshold):
    dist_size = np.size(average_movingness, 0)

    # movingness = value > threshold
    movingness = np.zeros((dist_size, 1))
    for i in range(0, dist_size):
        movingness[i, 0] = 0 if average_movingness[i, 0] < threshold else 1

    return movingness