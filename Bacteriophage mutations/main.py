import dataCollection
import Extraction
import writer
import os


if __name__ == '__main__':
    # make output directories
    os.mkdir('aa_mutations')
    os.mkdir('codon_frequency')
    os.mkdir('nuc_mutations')
    os.mkdir('nucleotides')
    os.mkdir('phages')
    # call collection functions
    collect = dataCollection.DataCollection()
    phages, fastas = collect.cluster_phages()
    # call extraction functions
    extracting = Extraction.Extraction()
    extracting.codon_frequency()
    extracting.nucleotide_mutations()
    extracting.aa_mutations()
    # call writing functions
    writer.mutation_to_csv('Nucleotide')
    writer.mutation_to_csv('Amino Acid')
    writer.write_CAI()
    writer.write_gen_info()
    # clean
    dataCollection.clean()
