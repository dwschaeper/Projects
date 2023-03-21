import sys
import igraph as ig
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu
import warnings

warnings.filterwarnings('ignore')


def graph_detail(graph, plot):
    # create dataframe, and store results in it
    details = pd.DataFrame(columns=['Node_name', 'Degree', 'Clustering_coefficient'])
    nodes = graph.vs()['name']
    details['Node_name'] = nodes
    details['Degree'] = graph.degree()
    details['Clustering_coefficient'] = graph.transitivity_local_undirected()
    # output details in case anyone would like to see it
    details.to_csv('graph_details.csv')

    # compute average clustering_coeffiecnt
    avg_clustering_coefficient = details.mean()['Clustering_coefficient']
    print('Average clustering coefficient is: ', avg_clustering_coefficient)

    # get the degree distribution
    degree_distribution = details.value_counts(subset=['Degree']).sort_index().to_dict()
    x = []
    y = []
    for entry in degree_distribution:
        x.append(int(str(entry).split(',')[0][1:]))
        y.append(degree_distribution[entry])

    # double log transform
    x = np.log10(x)
    y = np.log10(y)
    # get the line of best fit
    m, b = np.polyfit(x, y, 1)
    if 2 < abs(m) < 3:
        print('This is scale free ({})'.format(abs(m)))
    else:
        print('This is not scale free ({})'.format(abs(m)))
    # optionally plot the figure
    if plot:
        plt.scatter(x, y)
        plt.plot(x, m * np.array(x) + b, color='red')
        plt.xlabel('Log(Degree)')
        plt.ylabel('Log(Number of nodes)')
        plt.show()


def short_path_length(graph, handle):
    # initialize variables
    proteins = []
    shortest_paths = []
    not_reached = []

    # open the file to get a list of proteins
    with open(handle, 'r') as file:
        for line in file.readlines():
            proteins.append(line.strip())
    # get all pairs in the list
    pairs1 = [(a, b) for index, a in enumerate(proteins) for b in proteins[index + 1:]]
    for pair in pairs1:
        # calculate the shortest path. If a protein is not in the graph, catch the error and continue
        try:
            shortest = graph.get_shortest_paths(pair[0], to=pair[1])[0]
            if len(shortest) > 0:
                shortest_paths.append(min(shortest))
            else:
                not_reached.append(pair)
        except ValueError as e:
            continue

    return shortest_paths, not_reached


if __name__ == '__main__':
    try:
        if sys.argv[1] == 'plot':
            plot = True
        else:
            print('Invalid argument')
            exit(1)
    except IndexError:
        plot = False
    human_graph_file = 'Human-PPI.txt'
    # read as dataframe and import as graph
    graph_frame = pd.read_csv(human_graph_file, sep='\t', skiprows=1, names=['OFFICIAL_SYMBOL_A', 'OFFICIAL_SYMBOL_B'])
    graph = ig.Graph.DataFrame(graph_frame)
    graph_detail(graph, plot)
    short1, unreached1 = short_path_length(graph, 'protein-list1.txt')
    short2, unreached2 = short_path_length(graph, 'protein-list2.txt')

    # using Mann-Whitney U because Wilcoxon requires same length for x and y. This is the standard adjustment in R for
    # the Wilcoxon Sign Test when using unequal length x and y so it was used here as well.
    stat, p_value = mannwhitneyu(short1, short2)
    print('P-value from test: ', p_value)
