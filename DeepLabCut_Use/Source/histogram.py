import numpy as np
from matplotlib import pyplot


def display_histogram(dist, bin):
    # calcul min hist
    Min = min(dist)
    Max = max(dist)

    Dif = (Max - Min) / bin
    VectBin = np.arange(bin)

    Vect = np.zeros((1, bin))
    for x in range(0, np.size(VectBin, 0)):
        Vect[0, x] = np.multiply(VectBin[x], Dif)

    # remplace the number by the VectBin
    Vect2 = np.zeros((1, np.size(dist, 0)))
    for x in range(0, np.size(dist, 0) - 1):
        for y in range(0, np.size(VectBin, 0) - 1):
            if Vect[0, y] < dist[x, 0] <= Vect[0, y + 1]:
                Vect2[0, x] = VectBin[y]

    Test = np.unique(Vect2, return_index=False, return_inverse=False, return_counts=True, axis=None)
    Test1 = np.asarray(Test)

    pyplot.bar(Test1[0, :], Test1[1, :])
    pyplot.xlabel('bin')
    pyplot.ylabel('difference of pixel')
    pyplot.title('Histogram of the repartition of the differences of pixel')

    pyplot.show()