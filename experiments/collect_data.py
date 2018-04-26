import os.path
import pandas as pd

if __name__ == '__main__':

    data = []
    for jobindex in range(1, 32):
        data.append(pd.read_csv("linnea_generation{}.csv".format(jobindex)))

    result = pd.concat(data)
    result.to_csv("linnea_generation.csv")