# 2025 Scripts for boxplots from human/Influenza RNA-Seq data

A collection of input tables and scripts to produce various gene expression boxplots based on [RNAflow](https://github.com/hoelzer-lab/rnaflow) output for studying human host and influenza gene expression. 

As input for the scripts, we use 

* TPM values from RNAflow for expression levels (via `featureCounts` and a custom script for TPM calculation)
* adjusted p-values for significance between pairwise comparisons from RNAflow (via `DESeq2`)

**ATTENTION: If you dont find a gene by name, check https://www.ensembl.org/index.html for the name and the corresponding ENSG ID. The ID should be in the files and then you can add/change the name in the TPM _and_ pvalue file to plot it** 

Below we documented the scripts used for preparing the input tables and the scripts for plotting. We also provide the prepared input files in this repository for reproducibility. 

**Please note** that the scripts will not work out-of-the-box when cloning this repository. Paths need to be adjusted and dependencies (such as `ruby`, `python3` and `panda`) installed. However, the scripts provide information on how the data was wrangled and the plots produced. 

## Expression boxplots for human genes

Due to the experimental protocol specificities, the reads were mapped reverse complementary to the human reference genome and annotation (we controlled by looking at the mapping results in IGV). Thus, we configured `RNAflow` to count the reads reverse complementary after mapping (`--strand 2`).

```bash
# combine the single TPM value files from RNAflow and add gene names to the output table
ruby scripts/combine-tpms.rb # results in input-data/tpms-human.tsv

# generate a table with adjusted p-values for pairwise comparisons (mock-avian, mock-swine, mock-reass)
ruby scripts/combine-pvalues.rb # results in pvals-all-human-comparisons.tsv
cp pvals-all-human-comparisons.tsv pvals-virus-comparisons.tsv # remove the vs. mock comparisons bc we decided to only show the vs virus comparisons, vs mock is anyway always significant

# do the plotting. ATTENTION: adjust the plotting scripts themself for which gene IDs to plot and how!
conda create -n pandas -c conda-forge pandas && conda activate pandas
# plot single box plots per gene ID
python scripts/boxplot-tpm-adjp.py
# plot a grid view box plot for many gene IDs
python scripts/boxplot-tpm-adjp-summarized.py
```

## Expression boxplots for Influenza segments

For Influenza, we counted reads mapping to the `+` and `-` strand separately to investigate mRNA(+)and cRNA(+) as well as vRNA(-) expression separately. We also did some slight modifications to the TPM counting tables to remove very likely artifact read counts (a count of one read, for example, spoiling the visualization). We also renamed some segments to better fit the experimental setup. See our paper for details.  

### With `--strand 2` counting, reads mapping to (-) strand

Personal note: see for raw data `/Users/martin/projects/2024-04-16-manual-DESeq-Agustina/2025-04-08/results`

```bash
# Create a table of TPM normalized epxression values for the eight segments (similar to the TPM table above)
ruby scripts/combine-counts-segments.rb # results in TPMs/counts-tpm-segments.tsv
cp input-data/counts-tpm-segments.tsv input-data/counts-tpm-remove-HA-mock2-count1-segments.tsv # manually modified, bc count 1 only for mock and rep2 for HA gene...  + renamed N1 to NA and M1 to MP

# generate a table with adjusted p-values for all pairwise comparisons (mock-avian, mock-swine, mock-reass, avian-swine, ...) 
ruby scripts/combine-pvalues-segments.rb
awk '{print $1"\t"$2"\t"$5"\t"$7}' pvals-segments.tsv > pvals-without-mock-segments.tsv # + renamed N1 to NA and M1 to MP

conda activate pandas
# change the input files and output in the script! 
python scripts/boxplot-tpm-adjp.py
```

### With `--strand 1` counting, reads mapping to (-) strand

Personal note: calculated the `DESeq2` results again similar to `/Users/martin/projects/2024-04-16-manual-DESeq-Agustina/2025-04-08/` New results are here: `/Users/martin/projects/2024-04-16-manual-DESeq-Agustina/2025-05-06-deseq2-for-strand1`

```bash
# Create a table of TPM normalized epxression values for the eight segments 
ruby scripts/combine-counts-segments-strand1.rb # results in TPMs/counts-tpm-segments-strand1.tsv
cp counts-tpm-segments-strand1.tsv counts-tpm-remove-NP_mock1_and_PB1_mock3-count1-segments-strand1.tsv # bc count 1 only for mock in rep1 and 3 for NP and PB1 genes...  + renamed N1 to NA

# generate a table with adjusted p-values for all pairwise comparisons (mock-avian, mock-swine, mock-reass, avian-swine, ...) 
ruby scripts/combine-pvalues-segments-strand1.rb # results in pvals-segments-strand1.tsv
cd input-data
awk '{print $1"\t"$2"\t"$5"\t"$7}' pvals-segments-strand1.tsv > pvals-without-mock-segments-strand1.tsv # + renamed N1 to NA and M1 to MP

conda activate pandas
# change the input files and output in the script:
python scripts/boxplot-tpm-adjp.py
```

### With `--strand 1` counting, reads mapping to (-) strand + only run DESeq2 on the segments w/o including the human strand 2 data for normalization!

Personal note: New `DESeq2` results are here: `/Users/martin/projects/2024-04-16-manual-DESeq-Agustina/2025-05-06-deseq2-for-strand1/results-strand1-only-segments-input`

```bash
# KEEP THE SAME TPM VALUES: but really the quesiton if they make so much sense... they are calculated by RNAflow together with counting --strand 1 reads for human genes which is not the correct counting for human (probably counted very few reads for human) 
# Create a table of TPM normalized epxression values for the eight segments 
#ruby scripts/combine-counts-segments-strand1.rb # results in TPMs/counts-tpm-segments-strand1.tsv
#cp counts-tpm-segments-strand1.tsv counts-tpm-remove-NP_mock1_and_PB1_mock3-count1-segments-strand1.tsv # bc count 1 only for mock in rep1 and 3 for NP and PB1 genes...  + renamed N1 to NA

# generate a table with adjusted p-values for all pairwise comparisons (mock-avian, mock-swine, mock-reass, avian-swine, ...) 
ruby scripts/combine-pvalues-segments-strand1-onlySegments.rb # results in pvals-segments-strand1-onlySegments.tsv
cd input-data
awk '{print $1"\t"$3"\t"$5"\t"$7}' pvals-segments-strand1-onlySegments.tsv > pvals-without-mock-segments-strand1-onlySegments.tsv # + renamed N1 to NA and M1 to MP

conda activate pandas
# change the input files and output in the script:
python scripts/boxplot-tpm-adjp.py
```

