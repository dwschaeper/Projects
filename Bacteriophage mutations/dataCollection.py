import requests
import json
import time
from Bio import Entrez, SeqIO
import shutil
import os
import urllib


class DataCollection:
    """=================================================================================================================
    Collects all the relevant information from PhagesDB (the phages, coding regions, accession number) and the
    nucleotide sequence from NCBI via Entrez.
    ================================================================================================================="""

    def __init__(self):
        """=============================================================================================================
        Initialize the class object with class variables.

        ============================================================================================================="""
        self.phagesdb = 'https://phagesdb.org/api/genesbyphage/'
        self.cluster = input('Enter the cluster: ')
        self.phage = ''
        self.name_json = ''
        self.name_fasta = ''

    def cluster_phages(self):
        """=============================================================================================================
        This function reads the file of phages and determines which phages are part of the cluster in question. If it
        is, it triggers the rest of the class functions.
        :return: True
        ============================================================================================================="""
        print('Collecting bacteriophage information...')
        # initialize phages list
        phagesdb_files = []
        fasta_files = []
        # attempt to open file
        try:
            file = open('PhagesDB_Data.txt', 'r')
        except (IOError, OSError):
            print('Error opening file')
            exit(1)

        # check each phage to see if it is in the proper cluster
        for line in file:
            info = line.split('\t')
            if info[1] == self.cluster:
                self.phage = info[0]
                print('\t', self.phage, '...')
                self.phage_information()
                self.nucleotides()
                phagesdb_files.append(self.name_json)
                fasta_files.append(self.name_fasta)
                # wait so as to not to abuse webservers
                time.sleep(.2)

        file.close()

        return phagesdb_files, fasta_files

    def phage_information(self):
        """=============================================================================================================
        This function gets the phage gene information and accession number.
        :return: True
        ============================================================================================================="""
        # send a get request to PhagesDB for the phage
        response = requests.get(self.phagesdb + self.phage).json()['results']
        # open/create a file for the phage and store the response
        self.name_json = 'phages/' + self.phage + '.json'
        with open(self.name_json, 'w') as file:
            json.dump(response, file)

        return True

    def nucleotides(self):
        """=============================================================================================================
        This function uses the accession number from the json output from PhagesDB and collects the nucleotide sequence.
        :return: True
        ============================================================================================================="""
        # initialize entrez email
        Entrez.email = 'daschaep@iu.edu'

        # initialize variables
        self.name_fasta = 'nucleotides/' + self.phage + '.fasta'
        # get the accession number
        with open(self.name_json, 'r') as file:
            data = json.load(file)
        accession = data[0]['PhageID']['Accession']
        if accession == '':
            os.remove(self.name_json)
        else:  # query NCBI
            try:
                with Entrez.efetch(db='nucleotide', id=accession, rettype='fasta', retmode='text', ) as ncbi_result:
                    seq = SeqIO.read(ncbi_result, 'fasta')
                    SeqIO.write(seq, self.name_fasta, 'fasta')
            except(urllib.error.HTTPError):
                os.remove(self.name_json)
        return True


def clean():
    """=============================================================================================================
    This function removes the temp json files.
    :return: True
    ============================================================================================================="""
    print('Cleaning temp files...')
    shutil.rmtree('aa_mutations')
    shutil.rmtree('codon_frequency')
    shutil.rmtree('nuc_mutations')
    shutil.rmtree('nucleotides')
    shutil.rmtree('phages')
    return True
