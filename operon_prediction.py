import pandas as pd
import warnings

warnings.filterwarnings('ignore')


class OperonPrediction:
    """=================================================================================================================
    This class predicts operons based upon a sorted ptt file.
    ================================================================================================================="""

    def __init__(self, file):
        """=============================================================================================================
        Initializer of the class.
        :param file: The ptt file to predict operons from
        ============================================================================================================="""
        if file[-4:] == '.ptt':
            self.data_frame = pd.read_csv(file, sep='\t', skiprows=2)
        elif file[-4:] == '.gff':
            self.data_frame = pd.read_csv(file, sep='\t', skiprows=1,
                                          names=['Contig', 'Version', 'Type', 'Start', 'Stop', '.', 'Strand', 'Zero',
                                                 'Details'])
        self.data_frame['Operon'] = ''

    def predict_ptt(self):
        """=============================================================================================================
        This function uses the ptt file to predict the operons.
        :return: The data frame with operons added
        ============================================================================================================="""
        # initialize variables
        prior = self.data_frame.iloc[0]
        length = int(self.data_frame.shape[0])
        count = 1
        self.data_frame['Operon'][0] = count
        # iterate through the dataframe
        for index in range(1, length):
            row = self.data_frame.iloc[index]

            # collect the strands
            row_strand = str(row['Strand'])
            prior_strand = str(prior['Strand'])

            # collect the start of the new and stop of the old gene
            row_start = int(str(row['Location']).split('..')[0])
            prior_end = int(str(prior['Location']).split('..')[1])

            if row_strand != prior_strand:
                count += 1
                self.data_frame['Operon'][index] = count
                prior = row
            elif row_start - prior_end >= 50:
                count += 1
                self.data_frame['Operon'][index] = count
                prior = row
            else:
                self.data_frame['Operon'][index] = count
                prior = row

        return self.data_frame

    def predict_gff(self):
        """=============================================================================================================
        This function uses a gff file to predict operons.
        :return: The dataframe with operons added
        ============================================================================================================="""
        # sort the dataframe
        sorted_frame = self.data_frame.sort_values(by=['Contig', 'Start'], ignore_index=True)

        # initialize variables
        prior = sorted_frame.iloc[0]
        length = int(sorted_frame.shape[0])
        count = 1
        sorted_frame['Operon'][0] = count

        for index in range(1, length):
            row = self.data_frame.iloc[index]

            # collet the contig names
            row_contig = row['Contig']
            prior_contig = prior['Contig']

            # collect the strands
            row_strand = str(row['Strand'])
            prior_strand = str(prior['Strand'])

            # collect the start of the new and stop of the old contig
            row_start = int(row['Start'])
            prior_end = int(prior['Stop'])

            if row_contig != prior_contig:
                count += 1
                sorted_frame['Operon'][index] = count
                prior = row
            elif row_strand != prior_strand:
                count += 1
                sorted_frame['Operon'][index] = count
                prior = row
            elif row_start - prior_end >= 50:
                count += 1
                sorted_frame['Operon'][index] = count
                prior = row
            else:
                sorted_frame['Operon'][index] = count
                prior = row

        return sorted_frame


# the main program
if __name__ == '__main__':
    # predict operons for the .ptt files
    for filename in ['B_subtilis_168.ptt', 'E_coli_K12_MG1655.ptt', 'Halobacterium_NRC1.ptt',
                     'Synechocystis_PCC6803_uid159873.ptt']:
        prediction = OperonPrediction(filename)
        data_frame = prediction.predict_ptt()
        handle = filename[0:-4] + 'Operons.txt'
        data_frame.to_string(buf=handle, index=False, max_rows=data_frame.shape[0], max_cols=data_frame.shape[1],
                             justify='left')

    # predict the operons for the .gff
    gff = OperonPrediction('2088090036.gff')
    data_frame = gff.predict_gff()
    handle = '2088090036Operons.txt'
    data_frame.to_string(buf=handle, index=False, max_rows=data_frame.shape[0], max_cols=data_frame.shape[1],
                         justify='left')
