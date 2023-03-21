import numpy as np
from math import inf
import pandas as pd


class MotifFindandSearch:
    def __init__(self):
        columns = [i for i in range(2, 20)]
        self.counts = np.loadtxt('argR-counts-matrix.txt', delimiter='\t', usecols=columns)
        column_names = {'ID': '', 'Seq': '', 'Score': ''}
        self.top_scores = pd.DataFrame(columns=column_names)

    def log_odds_matrix(self):
        # add the pseudocount
        psuedo = self.counts + 1
        # compute the total number of samples i.e. sum down a column and translate count matrix to frequency matrix
        total = np.sum(psuedo, axis=0)[0]
        frequency = psuedo / total
        # compute the log-odds of each position
        self.logodds = np.log10(frequency / .25)

    def prediction(self):
        # open the sequences file
        sequences = open('E_coli_K12_MG1655.400_50 (3)', 'r')
        # iterate through the sequences
        for line in sequences.readlines():
            # initialize score max and access the id and seq
            max = -inf
            info = line.split(' \\ ')
            seq = info[1][0:-3]
            # itertate through the sequence to get all windows of 18 length
            for i in range(0, len(seq) - (self.logodds.shape[1]-1)):
                portion = seq[i: i + 18]
                score = 0
                # based on the nucleotide and position, sum the score
                for j in range(0, len(portion)):
                    if portion[j] == 'a':
                        score = score + self.logodds[0:1,j:j+1]
                    elif portion[j] == 'c':
                        score = score + self.logodds[1:2, j:j + 1]
                    elif portion[j] == 'g':
                        score = score + self.logodds[2:3, j:j + 1]
                    elif portion[j] == 't':
                        score = score + self.logodds[3:4, j:j + 1]
                # update max score
                if score > max:
                    max = score
                    id_max = info[0]
                    seq_max = portion
            # append the max score to the bottom of the dataframe
            self.top_scores.loc[len(self.top_scores)] = id_max, seq_max, float(max)
        sequences.close()

        # sort the dataframe by score
        self.top_scores.sort_values(by=['Score'], inplace=True, ascending=False)
        self.top_scores.head(30).to_csv('top_hits.csv', index=False)


if __name__ == '__main__':
    PWM = MotifFindandSearch()
    PWM.log_odds_matrix()
    PWM.prediction()
