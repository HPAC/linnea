import os.path
import pandas as pd

if __name__ == '__main__':

    data = []
    for jobindex in range(31):
        data.append(pd.read_pickle("linnea_generation{}.pkl".format(job_index)))

    result = pd.concat(data)
    result.to_csv("linnea_generation.csv")