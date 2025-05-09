import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def plot_tpm_boxplot(tpm_file, gene_id):
    # Load the TPM table
    df = pd.read_csv(tpm_file, sep="\t")

    # Check if gene exists
    if gene_id not in df['ID'].values:
        print(f"Gene ID '{gene_id}' not found in the TPM table.")
        return

    # Extract the row for the selected gene
    gene_data = df[df['ID'] == gene_id]

    # Extract TPM values for each condition (assuming columns follow naming like mock_1, avian_1, etc.)
    data = {
        'Mock': gene_data[['mock-rep1', 'mock-rep2', 'mock-rep3']].values.flatten(),
        'Avian': gene_data[['avian-rep1', 'avian-rep2', 'avian-rep3']].values.flatten(),
        'Swine': gene_data[['swine-rep1', 'swine-rep2', 'swine-rep3']].values.flatten(),
        'Reassortant': gene_data[['reass-rep1', 'reass-rep2', 'reass-rep3']].values.flatten()
    }

    # Convert to long-form DataFrame for seaborn
    plot_df = pd.DataFrame([
        {'Condition': cond, 'TPM': tpm}
        for cond, values in data.items()
        for tpm in values
    ])

    # Create the boxplot
    plt.figure(figsize=(4, 4)) #width, heigth
    sns.boxplot(x='Condition', y='TPM', data=plot_df, order=['Mock', 'Avian', 'Swine', 'Reassortant'], palette='Set2')
    sns.stripplot(x='Condition', y='TPM', data=plot_df, color='black', size=5, jitter=True, order=['Mock', 'Avian', 'Swine', 'Reassortant'])
    #plt.yscale('log')  # Log scale as seen in your example
    plt.ylabel('Expression level (TPM)')
    plt.title(f'Expression of {gene_id}')
    plt.tight_layout()
    plt.savefig(f'boxplot.{gene_id}.pdf')
    plt.show()

# Example usage:
# plot_tpm_boxplot('tpms.tsv', 'ENSG00000225855')

tpm_file = '/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/TPMs/tpms.tsv'
gene_id = 'ENSG00000225855'

plot_tpm_boxplot(tpm_file, gene_id)
