import os
import json
from Bio import SeqIO

def mutation_to_csv(unit_type):
    """=================================================================================================================
    Writes the input type of mutations to their specific .csv output files.
    :param unit_type: The type of mutation, either Nucleotide or Amino Acid
    :return:
    ================================================================================================================="""
    # check the type to write
    print('Writing', unit_type, 'Mutations...')
    if unit_type == 'Nucleotide':
        # get a list of all the files for amino acid information
        files = os.listdir('nuc_mutations')
        # initialize nuc specific variables
        length = 6
        handle_prefix = 'nuc_mutations/'
        # open the csv's to use as the inputs into the database
        base_output = open('Nucleotide_Mutations.csv', 'a')
        protein_output = open('Nucleotide_by_protein.csv', 'a')

        # check if the csv's are empty and if so place a header at the start
        if os.stat('Nucleotide_Mutations.csv').st_size == 0:
            header = 'Bacteriophage,A,C,G,T,insert,deletions\n'
            base_output.write(header)
        if os.stat('Nucleotide_by_protein.csv').st_size == 0:
            header = 'Bacteriophage,Protein,A,C,G,T,insert,deletions\n'
            protein_output.write(header)
    elif unit_type == 'Amino Acid':
        # get a list of all the files for amino acid information
        files = os.listdir('aa_mutations')
        # initialize aa specific variables
        length = 22
        handle_prefix = 'aa_mutations/'
        # open the csv's to use as the inputs into the database
        base_output = open('AA_Mutations.csv', 'a')
        protein_output = open('AA_by_protein.csv', 'a')

        # check if the csv's are empty and if so place a header at the start
        if os.stat('AA_Mutations.csv').st_size == 0:
            header = 'Bacteriophage,A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V, insert,deletions\n'
            base_output.write(header)
        if os.stat('AA_by_protein.csv').st_size == 0:
            header = 'Bacteriophage,Protein,A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V, insert,deletions\n'
            protein_output.write(header)
    else:
        print('Enter a valid type, either Nucleotide or Amino Acid')
        exit(1)

    # iterate through the files
    for file in files:
        # read each file
        handle = handle_prefix + file
        aa_counts = open(handle, 'r')
        insert = str(file.split(' ')[0]) + ','
        # place the counts in a comma separated format and write to the file for mutations
        for i in range(0, length):
            acid_count = aa_counts.readline().strip().split(' ')[1]
            insert = insert + str(acid_count) + ','
        insert = insert[0:-1] + '\n'
        base_output.write(insert)

        # place counts by protein in the csv for proteins
        for j in range(length, sum(1 for line in open(handle))):
            insert = str(file.split(' ')[0]) + ','
            line = aa_counts.readline().strip().split('{')
            insert = insert + line[0][0:-1] + ','
            for acid in str(line[1]).split(','):
                insert = insert + acid.split(':')[1].replace('}', '') + ','
            insert = insert + '\n'
            protein_output.write(insert)

        aa_counts.close()

    base_output.close()
    protein_output.close()


def write_CAI():
    """"=================================================================================================================
    Writes the information from the CAI files to CAI.csv.

    :return: True
    ================================================================================================================="""
    # get a list of all possible tri-nucleotides
    print('Writing CAI file')
    combos = []
    for base1 in ['T', 'C', 'A', 'G']:
        for base2 in ['T', 'C', 'A', 'G']:
            for base3 in ['T', 'C', 'A', 'G']:
                combo = base1 + base2 + base3
                combos.append(combo)
    # list the files
    files = os.listdir('codon_frequency')
    cai_csv = open('CAI.csv', 'a')
    # make the header if never used before
    if os.stat('CAI.csv').st_size == 0:
        header = 'Bacteriophage,'
        for base1 in ['T', 'C', 'A', 'G']:
            for base2 in ['T', 'C', 'A', 'G']:
                for base3 in ['T', 'C', 'A', 'G']:
                    header = header + base1 + base2 + base3 + ','
        header = header + '\n'
        cai_csv.write(header)

    # iterate through the files
    wrote = False
    for file in files:
        writing = file.split('.')[0] + ','
        handle = 'codon_frequency/' + file
        cai = open(handle, 'r')
        information = cai.readline()[1:-1].split(',')
        for tri_nuc in combos:
            for entry in information:
                info = entry.split(':')
                index = info[0][1:-1]
                value = info[1][1:]
                # tag if a null space needs to be used
                if index.replace("'", '') == tri_nuc:
                    wrote = True
                    writing = writing + value + ','
            if not wrote:
                writing = writing + ','
            wrote = False
        writing = writing + '\n'
        cai_csv.write(writing)
        cai.close()

    cai_csv.close()

    return True

def write_gen_info():
    """=================================================================================================================
    Writes the general information from the phagesDB response and length information from the Entrez fasta response to
    gen_info.csv.
    :return: True
    ================================================================================================================="""
    print('Writing general information')
    # read files for gen info
    gen_files = os.listdir('phages')
    # iterate through the files
    for file in gen_files:
        # read the json file
        json_handle = 'phages/' + file
        with open(json_handle, 'r') as gen_info:
            data = json.load(gen_info)
        # collect the information to write
        name = data[0]['PhageID']['PhageID']
        host = data[0]['PhageID']['HostStrain']
        accession = data[0]['PhageID']['Accession']
        cluster = data[0]['PhageID']['Cluster']
        fasta_handle = 'nucleotides/' + file[0:-5] + '.fasta'
        # get length info from fasta files
        fasta = SeqIO.read(fasta_handle, 'fasta')
        nuc_length = len(fasta.seq)
        aa_length = nuc_length / 3
        writing = name + ',' + host +',' + accession + ',' + cluster + ',' + str(nuc_length) + ',' + str(aa_length) + '\n'
        # check if output is empty, fill header if so, and write data to csv
        output = open('gen_info.csv', 'a')
        if os.stat('gen_info.csv').st_size == 0:
            header = 'Bacteriophage,Host,Accession,Cluster,NucleotideLength,AminoAcidLength\n'
            output.write(header)
        output.write(writing)

    output.close()

    return True