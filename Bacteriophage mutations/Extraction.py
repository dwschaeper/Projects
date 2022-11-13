import os
import json
from Bio import SeqIO
from Bio import pairwise2
from Bio import Align
from Bio.Align import substitution_matrices
from CAI import RSCU


class Extraction:
    """=================================================================================================================
    This class extracts the information for the database from the information collected.
    ================================================================================================================="""

    def __init__(self):
        """=============================================================================================================
        Class initializer
        ============================================================================================================="""
        self.fastas = os.listdir('nucleotides')
        self.phages = os.listdir('phages')
        self.known_genes = []

    def codon_frequency(self):
        """=============================================================================================================
        This function gets the codon frequency per bacteriophage
        :return: True
        ============================================================================================================="""
        # initialize holder
        all_regions = []

        # iterate through all of the phages collected in the collection step
        for index in range(len(self.phages)):
            phage_name = 'phages/' + self.phages[index]
            nuc_name = 'nucleotides/' + self.fastas[index]
            with open(phage_name, 'r') as data:
                self.file_genes = json.load(data)
            self.nucleotides = str(SeqIO.read(nuc_name, 'fasta').seq)

            # iterate through all of the genes in the phage and collect the coding region
            for gene in self.file_genes:
                start = int(gene['Start'])
                stop = int(gene['Stop'])
                gene_region = self.nucleotides[start:stop]
                # append the region to a list of all the regions if it is divisible by 3
                if len(gene_region) % 3 == 0:
                    all_regions.append(gene_region)
            # calculate the RSCU of all of the genes in the phage
            try:
                print('RSCU for', self.phages[index])
                rscu = RSCU(all_regions)
                rscu_name = 'codon_frequency/' + str(gene['PhageID']['PhageID']) + '.txt'
                # write the RSCU values for the phage in a file
                with open(rscu_name, 'w') as output:
                    output.write(str(rscu))

                # call the next function to continue
                self.functionally_known_coding_regions()
            except(ValueError):
                print('exception')
                # call the next function to continue
                # self.functionally_known_coding_regions()

        return True

    def functionally_known_coding_regions(self):
        """=============================================================================================================
        Determines information about nucleotide mutations for known genes
        :return: True
        ============================================================================================================="""
        phage_genes = []
        # collect gene names
        print('\tgene regions')
        for gene in self.file_genes:
            gene_name = str(gene['Notes']).split("'")
            # only if the gene is functionally known use it for mutations
            if gene_name[1] != '':
                start = int(gene['Start'])
                stop = int(gene['Stop'])
                gene_region = self.nucleotides[start:stop]
                regions_dictionary = {'name': gene['PhageID']['PhageID'],
                                      'cluster': gene['PhageID']['Cluster'],
                                      'start': start,
                                      'stop': stop,
                                      'function': gene_name[1],
                                      'nucleotides': gene_region,
                                      'amino acids': gene['translation']}
                phage_genes.append(regions_dictionary)

        self.known_genes.append(phage_genes)

        return True

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    def nucleotide_mutations(self):
        """=============================================================================================================
        This function determines information about nucleotide mutations across functionally known genes in the phage.
        :return: True
        ============================================================================================================="""
        counter = 1
        already_done = {}
        for target_phage in self.known_genes:
            # initialize variables for each phage
            name = target_phage[0]['name']
            print('Mutations for', name, '(', counter, len(self.known_genes), ')')
            handle = 'nuc_mutations/' + name + ' nucleotide mutations.txt'
            file = open(handle, 'w')
            by_function = {}
            mutations = {'A': 0,
                         'T': 0,
                         'G': 0,
                         'C': 0,
                         'insert': 0,
                         'deletion': 0}
            # select one gene at a time
            for row_target in target_phage:
                # select the phage to query against through all of the phages
                for query_phage in self.known_genes:
                    # iterate through the query phage to select individual gene regions
                    for row_query in query_phage:
                        # only align against genes that are in a different phage, as well as have the same function
                        if row_query['name'] != row_target['name'] and row_query['function'] == row_target[
                            'function'] and row_query['cluster'] == row_target['cluster']:
                            combo = row_query['name'], row_query['start'], row_query['stop'], row_target['name'], row_target['start'], row_target['stop']
                            # initialize counter, make sure to not overwrite information for genes with multiple copies
                            if row_target['function'] in by_function.keys():
                                function_mutations = {'A': by_function[row_target['function']]['A'],
                                                      'T': by_function[row_target['function']]['T'],
                                                      'G': by_function[row_target['function']]['G'],
                                                      'C': by_function[row_target['function']]['C'],
                                                      'insert': by_function[row_target['function']]['insert'],
                                                      'deletion': by_function[row_target['function']]['deletion']}
                            else:
                                function_mutations = {'A': 0,
                                                      'T': 0,
                                                      'G': 0,
                                                      'C': 0,
                                                      'insert': 0,
                                                      'deletion': 0}
                            # perform the alignment
                            if combo in already_done.keys():
                                print('this one done already')
                                mutations['deletion'] += already_done[combo]['deletion']
                                function_mutations['deletion'] += already_done[combo]['deletion']
                                mutations['insert'] += already_done[combo]['insert']
                                function_mutations['insert'] += already_done[combo]['insert']
                                mutations['A'] += already_done[combo]['A']
                                function_mutations['A'] += already_done[combo]['A']
                                mutations['T'] += already_done[combo]['T']
                                function_mutations['T'] += already_done[combo]['T']
                                mutations['G'] += already_done[combo]['G']
                                function_mutations['G'] += already_done[combo]['G']
                                mutations['C'] += already_done[combo]['C']
                                function_mutations['C'] += already_done[combo]['C']

                            else:
                                align = pairwise2.align.globalms(row_target['nucleotides'], row_query['nucleotides'], 2,
                                                                 -1,
                                                                 -1, -.5, one_alignment_only=True)
                                try:
                                    align_name = row_target['name'], row_target['start'], row_target['stop'], row_query[
                                        'name'], row_query['start'], row_query['stop']
                                    stralign = str(align[0])
                                    align_a = stralign[stralign.find('seqA=') + 5:stralign.find('seqB') - 2]
                                    align_b = stralign[stralign.find('seqB=') + 5:stralign.find('score') - 2]
                                    # collect the information
                                    for index in range(len(align_a)):
                                        if align_a[index] == '-':
                                            mutations['deletion'] += 1
                                            function_mutations['deletion'] += 1
                                        elif align_b[index] == '-':
                                            mutations['insert'] += 1
                                            function_mutations['insert'] += 1
                                        elif align_a[index] == align_b[index]:
                                            continue
                                        else:
                                            mutations[align_a[index]] += 1
                                            function_mutations[align_a[index]] += 1

                                    by_function[row_target['function']] = function_mutations
                                    already_done[align_name] = function_mutations
                                except(IndexError):
                                    print('Skipped one nuc mut')
            # write the results to a file for each bacteriophage
            for key, value in mutations.items():
                writing = str(key) + ' ' + str(value) + '\n'
                file.write(writing)
            for key, value in by_function.items():
                writing = str(key) + ' ' + str(value) + '\n'
                file.write(writing)
            file.close()
            counter += 1

        return True

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    def aa_mutations(self):
        """=============================================================================================================
        This function determines amino acid mutations.
        :return: True
        ============================================================================================================="""
        # initialize aligner object
        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -.5
        counter = 1
        already_done = {}

        # iterate through the phages collected
        for target_phage in self.known_genes:
            # inititalize variables
            name = target_phage[0]['name']
            print('Amino acid mutations for', name, '(', counter, '/', len(self.known_genes), ')')
            handle = 'aa_mutations/' + name + ' aa mutations.txt'
            output = open(handle, 'w')
            by_function = {}
            aa_mutations = {'A': 0,
                            'R': 0,
                            'N': 0,
                            'D': 0,
                            'C': 0,
                            'Q': 0,
                            'E': 0,
                            'G': 0,
                            'H': 0,
                            'I': 0,
                            'L': 0,
                            'K': 0,
                            'M': 0,
                            'F': 0,
                            'P': 0,
                            'S': 0,
                            'T': 0,
                            'W': 0,
                            'Y': 0,
                            'V': 0,
                            'insert': 0,
                            'deletion': 0}
            # select just one gene
            for row_target in target_phage:
                # select a phage to query against
                for query_phage in self.known_genes:
                    # select just one gene at a time
                    for row_query in query_phage:
                        # only work with regions that are from different phages of the same function
                        if row_query['name'] != row_target['name'] and row_query['function'] == row_target['function']:
                            combo = row_query['name'], row_query['start'], row_query['stop'], row_target['name'], \
                                    row_target['start'], row_target['stop']
                            if row_target['function'] in by_function.keys():
                                function_mutations = {'A': by_function[row_target['function']]['A'],
                                                      'R': by_function[row_target['function']]['R'],
                                                      'N': by_function[row_target['function']]['N'],
                                                      'D': by_function[row_target['function']]['D'],
                                                      'C': by_function[row_target['function']]['C'],
                                                      'Q': by_function[row_target['function']]['Q'],
                                                      'E': by_function[row_target['function']]['E'],
                                                      'G': by_function[row_target['function']]['G'],
                                                      'H': by_function[row_target['function']]['H'],
                                                      'I': by_function[row_target['function']]['I'],
                                                      'L': by_function[row_target['function']]['L'],
                                                      'K': by_function[row_target['function']]['K'],
                                                      'M': by_function[row_target['function']]['M'],
                                                      'F': by_function[row_target['function']]['F'],
                                                      'P': by_function[row_target['function']]['P'],
                                                      'S': by_function[row_target['function']]['S'],
                                                      'T': by_function[row_target['function']]['T'],
                                                      'W': by_function[row_target['function']]['W'],
                                                      'Y': by_function[row_target['function']]['Y'],
                                                      'V': by_function[row_target['function']]['V'],
                                                      'insert': by_function[row_target['function']]['insert'],
                                                      'deletion': by_function[row_target['function']]['deletion']}
                            else:
                                function_mutations = {'A': 0,
                                                      'R': 0,
                                                      'N': 0,
                                                      'D': 0,
                                                      'C': 0,
                                                      'Q': 0,
                                                      'E': 0,
                                                      'G': 0,
                                                      'H': 0,
                                                      'I': 0,
                                                      'L': 0,
                                                      'K': 0,
                                                      'M': 0,
                                                      'F': 0,
                                                      'P': 0,
                                                      'S': 0,
                                                      'T': 0,
                                                      'W': 0,
                                                      'Y': 0,
                                                      'V': 0,
                                                      'insert': 0,
                                                      'deletion': 0}
                            if combo in already_done.keys():
                                print('this one done already')
                                aa_mutations['deletion'] += already_done[combo]['deletion']
                                function_mutations['deletion'] += already_done[combo]['deletion']
                                aa_mutations['insert'] += already_done[combo]['insert']
                                function_mutations['insert'] += already_done[combo]['insert']
                                aa_mutations['A'] += already_done[combo]['A']
                                function_mutations['A'] += already_done[combo]['A']
                                aa_mutations['R'] += already_done[combo]['R']
                                function_mutations['R'] += already_done[combo]['R']
                                aa_mutations['N'] += already_done[combo]['N']
                                function_mutations['N'] += already_done[combo]['N']
                                aa_mutations['D'] += already_done[combo]['D']
                                function_mutations['D'] += already_done[combo]['D']
                                aa_mutations['C'] += already_done[combo]['C']
                                function_mutations['C'] += already_done[combo]['C']
                                aa_mutations['Q'] += already_done[combo]['Q']
                                function_mutations['Q'] += already_done[combo]['Q']
                                aa_mutations['H'] += already_done[combo]['H']
                                function_mutations['H'] += already_done[combo]['H']
                                aa_mutations['I'] += already_done[combo]['I']
                                function_mutations['I'] += already_done[combo]['I']
                                aa_mutations['L'] += already_done[combo]['L']
                                function_mutations['L'] += already_done[combo]['L']
                                aa_mutations['K'] += already_done[combo]['K']
                                function_mutations['K'] += already_done[combo]['K']
                                aa_mutations['M'] += already_done[combo]['M']
                                function_mutations['M'] += already_done[combo]['M']
                                aa_mutations['F'] += already_done[combo]['F']
                                function_mutations['F'] += already_done[combo]['F']
                                aa_mutations['P'] += already_done[combo]['P']
                                function_mutations['P'] += already_done[combo]['P']
                                aa_mutations['S'] += already_done[combo]['S']
                                function_mutations['S'] += already_done[combo]['S']
                                aa_mutations['T'] += already_done[combo]['T']
                                function_mutations['T'] += already_done[combo]['T']
                                aa_mutations['W'] += already_done[combo]['W']
                                function_mutations['W'] += already_done[combo]['W']
                                aa_mutations['Y'] += already_done[combo]['Y']
                                function_mutations['Y'] += already_done[combo]['Y']
                                aa_mutations['V'] += already_done[combo]['V']
                                function_mutations['V'] += already_done[combo]['V']

                            else:
                                # perform the alignment
                                align = aligner.align(row_target['amino acids'], row_query['amino acids'])

                                try:
                                    align_name = row_target['name'], row_target['start'], row_target['stop'], row_query[
                                        'name'], row_query['start'], row_query['stop']
                                    stralign = str(align[0])
                                    third = int(len(stralign) / 3)
                                    align_a = stralign[0:third]
                                    align_b = stralign[-third:]
                                    # collect the information
                                    for index in range(len(align_a)):
                                        if align_a[index] == '-':
                                            aa_mutations['deletion'] += 1
                                            function_mutations['deletion'] += 1
                                        elif align_b[index] == '-':
                                            aa_mutations['insert'] += 1
                                            function_mutations['insert'] += 1
                                        elif align_a[index] == align_b[index]:
                                            continue
                                        else:
                                            aa_mutations[align_a[index]] += 1
                                            function_mutations[align_a[index]] += 1

                                    by_function[row_target['function']] = function_mutations
                                    already_done[align_name] = function_mutations
                                except(IndexError):
                                    print('Skipped one aa mut')
            # write the results to a file for each bacteriophage
            for key, value in aa_mutations.items():
                writing = str(key) + ' ' + str(value) + '\n'
                output.write(writing)
            for key, value in by_function.items():
                writing = str(key) + ' ' + str(value) + '\n'
                output.write(writing)
            counter += 1
            output.close()
    # ------------------------------------------------------------------------------------------------------------------
