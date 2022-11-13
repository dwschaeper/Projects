import os
import pandas as pd
from numpy import log2
from sklearn.preprocessing import MinMaxScaler
from sklearn.feature_selection import mutual_info_classif
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.model_selection import GridSearchCV
from scipy.stats import ttest_ind


def common_mirna():
    print('Reading & Cleaning')
    # initialize variables
    types = ['Metastatic', 'Tumor']
    mirna_frame = pd.DataFrame()

    # get a dataframe for each and append to the total
    for cancer_type in types:
        cancers = os.listdir(cancer_type)
        for cancer in cancers:
            path = cancer_type + '/' + cancer
            for sample in os.listdir(path):
                # get how many lines to skip
                handle = path + '/' + sample
                with open(handle, 'r') as file:
                    counter = 0
                    for line in file:
                        if line[0:3] != 'hsa':
                            counter += 1
                        else:
                            # stop because no need to keep going only slows it down
                            break

                    this_frame = pd.read_csv(handle, skiprows=counter, sep='\t', index_col=0,
                                             names=['miRNA', 'value']).transpose()
                    # clean miRNAs
                    this_frame.drop([col for col in this_frame.columns if 'hsa' not in col], axis=1, inplace=True)

                    # take log2 (other data has the same transformation) of the data for lung cancer and match the
                    # naming convention
                    if cancer == 'Lung':
                        this_frame = log2(this_frame)
                        this_frame.columns = this_frame.columns.str.replace('_st', '')
                        this_frame.columns = this_frame.columns.str.replace('-star', '*')

                    # add cancer type
                    if cancer == 'Breast':
                        this_frame['Type'] = 0
                    elif cancer == 'Prostate':
                        this_frame['Type'] = 1
                    else:
                        this_frame['Type'] = 2
                    # add the binary label, 0 for normal tumor 1 for metastatic
                    if cancer_type == 'Metastatic':
                        this_frame['Result'] = 1
                    else:
                        this_frame['Result'] = 0

                    mirna_frame = pd.concat([mirna_frame, this_frame], axis=0, ignore_index=True)

        # remove all miRNA where not all samples have it
        mirna_frame.dropna(axis=1, inplace=True)
    mirna_frame.to_csv('dataframe.csv')
    return mirna_frame


def feature_selection(frame):
    print('Feature selection')
    # separate into result and feature
    results = frame['Result']
    types = frame['Type']
    features = frame.iloc[:, :-2]

    # scale the features
    scaler = MinMaxScaler()
    features[features.columns] = scaler.fit_transform(features[features.columns])
    # save the scaled values
    writing_features = features
    writing_features['Result'] = results
    writing_features.to_csv('heatmap.csv')
    # calculate information gain
    columns = frame.columns[0:-2]
    gain = pd.DataFrame(columns=['miRNA', 'IG'])
    information_gain = mutual_info_classif(features, results, random_state=1)
    for i in range(0, len(columns)):
        gain.loc[len(gain)] = [columns[i], information_gain[i]]

    # select top features
    gain.sort_values(by='IG', inplace=True, ascending=False, ignore_index=True)
    number = len(frame) // 10
    best = gain.head(number)
    # best = gain[gain['IG'] > 0.1]
    top_features = best['miRNA'].tolist()

    # subset frame to just the top features
    top_features.append('Result')
    top_features.append('Type')
    features['Result'] = results
    features['Type'] = types
    features = features[top_features]
    return features


def generate_svm(frame):
    print('Creating model')
    # split the data into test and train
    x = frame.iloc[:, :-2].to_numpy()
    y = frame['Result'].to_numpy()
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=.2, random_state=5)
    # create the model
    # TRY GRID SEARCH TO GET OPTIMAL C
    params = {'C': [0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.3, 1.7, 1.8, 2.0, 10, 15],
              'kernel': ['linear']}
    model = GridSearchCV(SVC(), param_grid=params, scoring='accuracy', verbose=3)
    model.fit(x_train, y_train)
    test_prediction = model.predict(x_test).tolist()

    # get accuracy for testing, and print it
    correct = 0
    for i in range(len(test_prediction)):
        if test_prediction[i] == y_test.tolist()[i]:
            correct += 1
    print(correct / len(test_prediction))

    # print AUC
    print(metrics.roc_auc_score(y_test, test_prediction))


def t_tests(frame):
    # split into two groups
    positive = frame[frame['Result'] == 1]
    negative = frame[frame['Result'] == 0]

    columns = positive.columns[:-2].tolist()

    # calculate all the p-values
    for column in columns:
        stat, p_value = ttest_ind(positive[column], negative[column])
        # print those that are significant
        if p_value < .05:
            print(p_value, column)
