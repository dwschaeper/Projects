import project as pr


if __name__ == '__main__':
    dataframe = pr.common_mirna()
    top_features = pr.feature_selection(dataframe)
    top_features.to_csv('features.csv')
    pr.generate_svm(top_features)
    pr.t_tests(top_features)

