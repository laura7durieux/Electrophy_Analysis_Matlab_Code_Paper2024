from enum import IntEnum

import numpy as np
import pandas as pd

# FORMAT 1
# class Column(IntEnum):
#     Nose = 0
#     Hat = 1
#     BackUp = 2
#     BackMiddle = 3
#     TailBase = 4
#     TailMiddle = 5
#     TailEnd = 6

# Format 2
class Column(IntEnum):
    Nose = 0
    Hat = 1
    LeftEar = 2
    RigthEar = 3
    LeftShoulder = 4
    RigthShoulder = 5
    UpperBack = 6
    MiddleBack = 7
    LeftHip = 8
    RigthHip = 9
    LowerBack = 10
    BaseTail = 11
    MiddleTail = 12
    # EndTail = 13

class Head(IntEnum):
    Nose = 0
    Hat = 1
    LeftEar = 2
    RigthEar = 3
    # LeftShoulder = 4
    # RigthShoulder = 5
    # UpperBack = 6
    # MiddleBack = 7
    # LeftHip = 8
    # RigthHip = 9
    # LowerBack = 10
    # BaseTail = 11
    # MiddleTail = 12
    # EndTail = 13

class DLCResult:
    def __init__(self, path):
        self.df = pd.read_csv(path)
        print(self.df)
        self.dataset = np.array([])

    def extract_coords(self, column: Column):
        first_column = 1 + (int(column) * 3);

        last_column_plus_one = first_column + 3

        raw_data = np.size(self.df, 0)  # size of the number of raw of the data
        dataset = self.df.iloc[1:raw_data, first_column:last_column_plus_one]  # keeping only the data position of the hat
        dataset.columns = ['x', 'y', 'likehood']  # changing the name of the columns
        dataset.reset_index()  # reseting the index of the data
        dataset = dataset.to_numpy()
        self.dataset = dataset
        return dataset


