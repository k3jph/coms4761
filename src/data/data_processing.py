import os
import sys
import pandas as pd

INPUT_ERROR_MSG = "Invalid input: please provide a network ID (3 for E. coli, 4 for S. cerevisiae"
NETWORK_ID = 3

def main():
    genes = pd.read_csv(f'raw/net{NETWORK_ID}_gene_ids.tsv', sep="\t")
    expressions = pd.read_csv(f'raw/net{NETWORK_ID}_expression_data.tsv', sep="\t")
    network = pd.read_csv(f'raw/DREAM5_NetworkInference_GoldStandard_Network{NETWORK_ID}.tsv', sep="\t", header=None)

    gene_dict = genes.set_index('#ID')['Name'].to_dict()

    # Network utilizes anonymized gene IDs for gene-gene relations
    # Map both TF and impacted genes to names
    network.columns = ['TF_ID', 'GENE_ID', 'IND']
    network['TF_ID'] = network['TF_ID'].map(gene_dict)
    network['GENE_ID'] = network['GENE_ID'].map(gene_dict)

    # Export in expected formats
    network.to_csv('Network.csv', index=False, header=False)
    expressions.to_csv('Expression.txt', index=False, header=False)
    genes['Name'].to_csv('Genes.csv', index=False, header=False)


if __name__ == '__main__':
    # Retrieve user input - if no network id number specified, default to 3 for E. coli
    args = sys.argv

    if len(args) > 1:
        try:
            NETWORK_ID = int(args[1])
        except:
            print(INPUT_ERROR_MSG)
            sys.exit()

    if not (NETWORK_ID == 3 or NETWORK_ID == 4):
        print(INPUT_ERROR_MSG)
    else:
        main()