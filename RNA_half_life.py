import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from math import ceil
from gprofiler.gprofiler import GProfiler


def calculate_half_life():
    # read in the data
    time_course_1 = np.genfromtxt('DecayTimecourse.txt', skip_header=2, usecols=[1, 2, 3, 4, 5, 6, 7, 8, 9],
                                  delimiter='\t',
                                  filling_values=-1)
    time_course_2 = np.genfromtxt('DecayTimecourse.txt', skip_header=2, usecols=[10, 11, 12, 13, 14, 15, 16, 17, 18],
                                  delimiter='\t', filling_values=-1)
    time_course_3 = np.genfromtxt('DecayTimecourse.txt', skip_header=2, usecols=[19, 20, 21, 22, 23, 24, 25, 26, 27],
                                  delimiter='\t', filling_values=-1)

    # generate an empty dataframe to fill with the results
    names = []
    file = open('DecayTimecourse.txt', 'r')
    file.readline()
    file.readline()
    for line in file.readlines():
        info = line.split('\t')
        names.append(info[0])
    file.close()
    half_life = pd.DataFrame(index=names)
    half_life['time_course_1'] = ""
    half_life['time_course_2'] = ""
    half_life['time_course_3'] = ""

    # perform linear regression to get the slope for time course 1
    time_stamps = [0, 5, 10, 15, 20, 30, 40, 50, 60]

    for row_num in range(len(time_course_1)):
        x = []
        y = []
        for index in range(len(time_course_1[row_num])):
            if time_course_1[row_num][index] == -1:
                continue
            else:
                x.append(time_stamps[index])
                y.append(time_course_1[row_num][index])
        if not x or len(x) == 1:
            half_life['time_course_1'][row_num] = None
        else:
            slope = LinearRegression().fit(np.array(x).reshape(-1, 1), np.array(y))
            half_life['time_course_1'][row_num] = float(slope.coef_)

    # perform linear regression to get the slope for time course 2
    for row_num in range(len(time_course_2)):
        x = []
        y = []
        for index in range(len(time_course_2[row_num])):
            if time_course_2[row_num][index] == -1:
                continue
            else:
                x.append(time_stamps[index])
                y.append(time_course_2[row_num][index])
        if not x or len(x) == 1:
            half_life['time_course_2'][row_num] = None
        else:
            slope = LinearRegression().fit(np.array(x).reshape(-1, 1), np.array(y))
            half_life['time_course_2'][row_num] = float(slope.coef_)

    # perform linear regression to get the slope for time course 3
    for row_num in range(len(time_course_3)):
        x = []
        y = []
        for index in range(len(time_course_3[row_num])):
            if time_course_3[row_num][index] == -1:
                continue
            else:
                x.append(time_stamps[index])
                y.append(time_course_3[row_num][index])
        if not x or len(x) == 1:
            half_life['time_course_3'][row_num] = None
        else:
            slope = LinearRegression().fit(np.array(x).reshape(-1, 1), np.array(y))
            half_life['time_course_3'][row_num] = float(slope.coef_)

    # create the average column and get final half life
    half_life['average'] = half_life.mean(axis=1)
    half_life['Final_Half_Life'] = np.log(2) / abs(half_life['average'])

    # sort and output
    sorted = half_life.sort_values(by=['Final_Half_Life'])
    sorted.to_csv('Half_Life.csv')

    # split into top and bottom 10% of those that have half lives. Some RNAs did not have 2 or more data points in at
    # least one of three time courses so half life for them was unable to be calculated
    RNAs_with_half_life = sorted.dropna(subset='Final_Half_Life')
    print(RNAs_with_half_life.shape)
    print(sorted.shape)
    ten_percent = ceil(len(RNAs_with_half_life) * .1)
    sorted_names = list(RNAs_with_half_life.index)
    top_ten = sorted_names[0:ten_percent]
    bottom_ten = sorted_names[-ten_percent:]

    # perform GO annotation
    top_gp = GProfiler(return_dataframe=True)
    top_results = top_gp.profile(organism='scerevisiae', query=top_ten)

    bottom_gp = GProfiler(return_dataframe=True)
    bottom_results = bottom_gp.profile(organism='scerevisiae', query=bottom_ten)

    # output result
    GO_result = pd.concat([top_results, bottom_results], ignore_index=True)
    GO_result.to_csv('GO_result.csv')


if __name__ == '__main__':
    calculate_half_life()
