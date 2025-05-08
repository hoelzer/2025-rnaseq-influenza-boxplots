import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def get_significance_stars(pval):
    if pval <= 0.001:
        return '***'
    elif pval <= 0.01:
        return '**'
    elif pval <= 0.05:
        return '*'
    else:
        return None  # Skip non-significant

def plot_tpm_boxplot(tpm_file, pval_file, gene_id):
    # Load TPM and p-value data
    df = pd.read_csv(tpm_file, sep="\t")
    pvals = pd.read_csv(pval_file, sep="\t")

    if gene_id not in df['ID'].values or gene_id not in pvals['ID'].values:
        print(f"Gene ID '{gene_id}' not found in one of the tables.")
        return

    gene_data = df[df['ID'] == gene_id]
    pval_data = pvals[pvals['ID'] == gene_id]

    # Get gene name or fallback to ID
    gene_name = gene_data['Name'].values[0] if 'Name' in gene_data.columns else gene_id
    if pd.isna(gene_name) or gene_name.strip() == "":
        gene_name = gene_id
    safe_name = gene_name.replace(" ", "_").replace("/", "_")
    filename_base = f'{safe_name}_{gene_id}'

    # TPM values per condition
    data = {
        'Mock': gene_data[['mock-rep1', 'mock-rep2', 'mock-rep3']].values.flatten(),
        'Avian': gene_data[['avian-rep1', 'avian-rep2', 'avian-rep3']].values.flatten(),
        'Swine': gene_data[['swine-rep1', 'swine-rep2', 'swine-rep3']].values.flatten(),
        'Reassortant': gene_data[['reass-rep1', 'reass-rep2', 'reass-rep3']].values.flatten()
    }

    # Prepare data for plotting
    plot_df = pd.DataFrame([
        {'Virus': cond, 'TPM': tpm}
        for cond, values in data.items()
        for tpm in values
    ])

    # Define color palette
    custom_palette = {
        'Mock': (1.0, 1.0, 1.0, 1.0),              # white
        'Reassortant': (0.6078, 0.7451, 0.9608, 1.0),  # blue
        'Avian': (0.9490, 0.6392, 0.6392, 1.0),       # red
        'Swine': (0.7137, 0.9804, 0.6784, 1.0)        # green
    }

    plt.figure(figsize=(4, 4))
    sns.boxplot(x='Virus', y='TPM', data=plot_df, palette=custom_palette, order=['Mock', 'Avian', 'Swine', 'Reassortant'])
    sns.stripplot(x='Virus', y='TPM', data=plot_df, color='black', size=5, jitter=True, order=['Mock', 'Avian', 'Swine', 'Reassortant'])

    # Check for log scale
    y_max = plot_df['TPM'].max()
    y_min = plot_df['TPM'].min()
    use_log_scale = y_max / max(y_min, 0.1) > 100

    if use_log_scale:
        plt.yscale('log')
        plt.ylabel('Expression level (TPM, log scale)')
    else:
        plt.ylabel('Expression level (TPM)')

    # Title: Gene name and ID
    plt.title(f'{gene_name} ({gene_id})')

    # Dynamic p-value annotation
    positions = {'Mock': 0, 'Avian': 1, 'Swine': 2, 'Reassortant': 3}
    cond_names = list(positions.keys())

    offset_counter = 0
    log_y_positions = []

    if use_log_scale:
        base_height = y_max * 1.5
    else:
        y_range = y_max - y_min
        base_height = y_max + 0.05 * y_range

    for col in pval_data.columns:
        if col == "ID":
            continue

        try:
            cond1_raw, cond2_raw = col.split("-")
            cond1 = cond1_raw.capitalize()
            cond2 = cond2_raw.capitalize()

            # Fix capitalization manually if needed
            cond1 = "Reassortant" if cond1.startswith("Reass") else cond1
            cond2 = "Reassortant" if cond2.startswith("Reass") else cond2

            if cond1 not in cond_names or cond2 not in cond_names:
                continue

        except Exception as e:
            print(f"Skipping {col}: {e}")
            continue

        pval = pval_data.iloc[0][col]
        stars = get_significance_stars(pval)
        if stars is None:
            continue  # Skip non-significant

        x1, x2 = positions[cond1], positions[cond2]

        if use_log_scale:
            if not log_y_positions:
                y = base_height
            else:
                y = log_y_positions[-1] * 2.2
            log_y_positions.append(y)
            bar_height = y * 0.15
        else:
            y = base_height + offset_counter * 0.12 * y_range
            bar_height = 0.02 * y_range

        text_label = f"p={pval:.2g} ({stars})"
        plt.plot([x1, x1, x2, x2], [y, y + bar_height, y + bar_height, y], lw=1.2, c='black')
        plt.text((x1 + x2) / 2, y + bar_height * 1.05, text_label, ha='center', va='bottom', fontsize=9)

        offset_counter += 1

    # Extend y-limits to fit all annotations
    if use_log_scale and log_y_positions:
        plt.ylim(y_min, log_y_positions[-1] * 2)
    elif not use_log_scale:
        current_ylim = plt.ylim()
        plt.ylim(current_ylim[0], current_ylim[1] + 0.05 * y_range)

    plt.tight_layout()

    # Save outputs
    plt.savefig(f'boxplot.{filename_base}.pdf')
    plt.savefig(f'boxplot.{filename_base}.png', dpi=300)
    plt.savefig(f'boxplot.{filename_base}.svg')
    # plt.show()  # Uncomment if you want interactive plots

#############
## Human genes
#tpm_file = '/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/TPMs/tpms.tsv'
#pval_file = '/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/TPMs/pvals.tsv' # this is old, bc we use now all virus comparisons except vs mock 
#pval_file = '/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/TPMs/pvals-all-human-comparisons.tsv' # this is old, bc we use now all virus comparisons except vs mock
#pval_file = '/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/TPMs/pvals-virus-comparisons.tsv'

#gene_id = 'ENSG00000225855'
#gene_id = 'ENSG00000107201' # DDX58
#gene_id = 'ENSG00000182393'	# IFNL1
#gene_id = "ENSG00000182782"

#plot_tpm_boxplot(tpm_file, pval_file, gene_id)

#TJP1    ENSG00000104067
#OCLN    ENSG00000197822
#CDH5    ENSG00000179776
#CDH4   ENSG00000179242
#CASP3   ENSG00000164305
#LDHA    ENSG00000134333
#LDHB    ENSG00000111716
#LDHC    ENSG00000166796
#LDHD    ENSG00000166816
#IFITM3  ENSG00000142089
#EIF2AK2 ENSG00000055332    (PKR) 

#gene_ids = ['ENSG00000104067','ENSG00000197822','ENSG00000179776','ENSG00000164305','ENSG00000134333','ENSG00000111716','ENSG00000166796','ENSG00000166816','ENSG00000142089','ENSG00000055332','ENSG00000179242']
#for gene_id in gene_ids:
#    plot_tpm_boxplot(tpm_file, pval_file, gene_id)


######### Genes for Supplement probably, downregulated in Reassortant

# ENSG00000197061 HIST1H4C
# ENSG00000184270 HIST2H2AB
# ENSG00000187840 EIF4EBP1
# ENSG00000029993	HMGB3
# ENSG00000168078	PBK
# ENSG00000205544	TMEM256
#gene_id = 'ENSG00000197061'
#gene_id = 'ENSG00000184270'
#gene_id = 'ENSG00000187840'
#gene_id = 'ENSG00000029993'
#gene_id = "ENSG00000168078"
#gene_id = "ENSG00000205544"

#plot_tpm_boxplot(tpm_file, pval_file, gene_id)

########### Genes for Fig3 F panel

#Then, I think Panel F should show some of the 15 genes top up-regulated in Panel C (other will be in Supp. material, like 
# those we want to point out from panel D and E). They are related to inmune response and is ok too. 
# I would use this ones: IFIT2, INFL1, CCL5, OASL, HSPA6, INFB1, MX1 
# (if are too many I would take out HSPA6 that is not well known as the others)

# IFIT2 ENSG00000119922
# IFNL1 ENSG00000182393 alias: IL29
# CCL5 ENSG00000271503
# OASL  ENSG00000135114
# HSPA6 ENSG00000173110
# INFB1 ENSG00000171855 alias: IFNB
# MX1   ENSG00000157601 (hello darkness my old friend ... https://pmc.ncbi.nlm.nih.gov/articles/PMC5512242/)

#gene_ids = ['ENSG00000119922', 'ENSG00000182393', 'ENSG00000271503', 'ENSG00000135114', 'ENSG00000173110', 'ENSG00000171855', 'ENSG00000157601']
#for gene_id in gene_ids:
#    plot_tpm_boxplot(tpm_file, pval_file, gene_id)




# Automatically extract all gene IDs from the pvals file
#pvals_df = pd.read_csv(pval_file, sep="\t")
#gene_ids = pvals_df['ID'].dropna().unique()

#for gene_id in gene_ids:
#    plot_tpm_boxplot(tpm_file, pval_file, gene_id)



##############
## Influenza Segments - STRAND 2
#tpm_file = '/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/TPMs/counts-tpm-remove-HA-mock2-count1-segments.tsv'
#pval_file = '/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/TPMs/pvals-without-mock-segments.tsv'

#gene_id = 'gene-PB2'
#gene_id = 'gene-HA'
#gene_id = 'gene-NS1'
#gene_id = 'gene-PA'

#gene_ids = ['gene-HA','gene-MP','gene-NA','gene-NP','gene-NS1','gene-PA','gene-PB1','gene-PB2']
#for gene_id in gene_ids:
#    plot_tpm_boxplot(tpm_file, pval_file, gene_id)



##############
## Influenza Segments - STRAND 1
tpm_file = '/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/TPMs/counts-tpm-remove-NP_mock1_and_PB1_mock3-count1-segments-strand1.tsv'
pval_file = '/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/TPMs/pvals-without-mock-segments-strand1.tsv'

#gene_id = 'gene-PB2'
#gene_id = 'gene-HA'
#gene_id = 'gene-NS1'
#gene_id = 'gene-PA'

gene_ids = ['gene-HA','gene-MP','gene-NA','gene-NP','gene-NS1','gene-PA','gene-PB1','gene-PB2']
for gene_id in gene_ids:
    plot_tpm_boxplot(tpm_file, pval_file, gene_id)




# COLOR CODES

# mock: white
# reassortant: blue, RGBA: 9bbef5ff
# avian: red, f2a3a3ff
# swine: green,  b6faadff