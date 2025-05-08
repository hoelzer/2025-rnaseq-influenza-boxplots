import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math

def get_significance_stars(pval):
    if pval <= 0.001:
        return '***'
    elif pval <= 0.01:
        return '**'
    elif pval <= 0.05:
        return '*'
    else:
        return None

def plot_tpm_boxplot(tpm_file, pval_file, gene_id, collect_data=False):
    df = pd.read_csv(tpm_file, sep="\t")
    pvals = pd.read_csv(pval_file, sep="\t")

    if gene_id not in df['ID'].values or gene_id not in pvals['ID'].values:
        print(f"Gene ID '{gene_id}' not found.")
        return

    gene_data = df[df['ID'] == gene_id]
    pval_data = pvals[pvals['ID'] == gene_id]

    gene_name = gene_data['Name'].values[0] if 'Name' in gene_data.columns else gene_id
    if pd.isna(gene_name) or gene_name.strip() == "":
        gene_name = gene_id
    safe_name = gene_name.replace(" ", "_").replace("/", "_")
    filename_base = f'{safe_name}_{gene_id}'

    data = {
        'Mock': gene_data[['mock-rep1', 'mock-rep2', 'mock-rep3']].values.flatten(),
        'Avian': gene_data[['avian-rep1', 'avian-rep2', 'avian-rep3']].values.flatten(),
        'Swine': gene_data[['swine-rep1', 'swine-rep2', 'swine-rep3']].values.flatten(),
        'Reassortant': gene_data[['reass-rep1', 'reass-rep2', 'reass-rep3']].values.flatten()
    }

    plot_df = pd.DataFrame([
        {'Virus': cond, 'TPM': tpm}
        for cond, values in data.items()
        for tpm in values
    ])

    if collect_data:
        plot_df['Gene'] = gene_name
        plot_df['GeneID'] = gene_id
        return plot_df, pval_data, gene_name

def plot_combined_boxplot(all_data_df, pvals_dict):
    custom_palette = {
        'Mock': (1.0, 1.0, 1.0, 1.0),
        'Reassortant': (0.6078, 0.7451, 0.9608, 1.0),
        'Avian': (0.9490, 0.6392, 0.6392, 1.0),
        'Swine': (0.7137, 0.9804, 0.6784, 1.0)
    }

    y_max = all_data_df['TPM'].max()
    y_min = all_data_df['TPM'].min()
    use_log_scale = y_max / max(y_min, 0.1) > 100

    all_data_df['Virus_Gene'] = all_data_df['Gene'] + "\n" + all_data_df['Virus']

    plt.figure(figsize=(max(8, len(all_data_df['Virus_Gene'].unique()) * 0.8), 6))
    ax = sns.boxplot(x='Virus_Gene', y='TPM', hue='Virus', data=all_data_df, palette=custom_palette, dodge=False)

    sns.stripplot(
        x='Virus_Gene', y='TPM', 
        data=all_data_df, 
        color='black', 
        size=4, jitter=True,
        dodge=False,
        ax=ax
    )

    if use_log_scale:
        plt.yscale('log')
        plt.ylabel('Expression level (TPM, log scale)')
    else:
        plt.ylabel('Expression level (TPM)')

    plt.xticks(rotation=45, ha='right')
    plt.title('Summary Boxplot of All Genes')
    plt.legend(title='Virus', bbox_to_anchor=(1.05, 1), loc='upper left')

    # p-value annotations
    offset_counter = 0
    bar_height_factor = 0.85
    virus_order = ['Mock', 'Avian', 'Swine', 'Reassortant']

    virus_gene_labels = all_data_df['Virus_Gene'].unique().tolist()
    label_to_pos = {label: i for i, label in enumerate(virus_gene_labels)}

    for gene, pval_data in pvals_dict.items():
        for col in pval_data.columns:
            if col == 'ID':
                continue
            try:
                cond1_raw, cond2_raw = col.split("-")
                cond1 = cond1_raw.capitalize()
                cond2 = cond2_raw.capitalize()
                if cond1.startswith('Reass'):
                    cond1 = 'Reassortant'
                if cond2.startswith('Reass'):
                    cond2 = 'Reassortant'
            except:
                continue

            if cond1 not in virus_order or cond2 not in virus_order:
                continue

            pval = pval_data.iloc[0][col]
            stars = get_significance_stars(pval)
            if stars is None:
                continue

            label1 = f"{gene}\n{cond1}"
            label2 = f"{gene}\n{cond2}"

            if label1 not in label_to_pos or label2 not in label_to_pos:
                continue

            x1 = label_to_pos[label1]
            x2 = label_to_pos[label2]

            # Calculate height
            y = y_max * (1.05 + offset_counter * bar_height_factor) if use_log_scale else y_max + (offset_counter * bar_height_factor * (y_max - y_min))

            ax.plot([x1, x1, x2, x2], [y, y + (y*0.05), y + (y*0.05), y], lw=1.2, c='black')
            ax.text((x1 + x2) / 2, y + (y*0.1), f"{stars}", ha='center', va='bottom', fontsize=9)

            offset_counter += 1

    plt.tight_layout()
    plt.savefig('combined_summary_boxplot_with_pvalues.pdf')
    plt.savefig('combined_summary_boxplot_with_pvalues.png', dpi=300)
    plt.savefig('combined_summary_boxplot_with_pvalues.svg')
    plt.show()

# === MAIN execution ===

tpm_file = '/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/TPMs/tpms.tsv'
pval_file = '/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/TPMs/pvals-virus-comparisons.tsv'

gene_ids = ['ENSG00000119922', 'ENSG00000182393', 'ENSG00000271503', 'ENSG00000135114', 'ENSG00000173110', 'ENSG00000171855', 'ENSG00000157601']

all_data = []
pvals_dict = {}

for gene_id in gene_ids:
    result = plot_tpm_boxplot(tpm_file, pval_file, gene_id, collect_data=True)
    if result:
        gene_data_df, pval_data, gene_name = result
        all_data.append(gene_data_df)
        pvals_dict[gene_name] = pval_data

all_data_df = pd.concat(all_data, ignore_index=True)
plot_combined_boxplot(all_data_df, pvals_dict)
