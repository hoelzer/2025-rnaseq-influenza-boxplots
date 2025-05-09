gene_ids = []
conditions = []
condition2tpm = {}
['mock', 'swine', 'avian', 'reass'].each do |condition|
    condition2tpm[condition] = {}
end

# get gene names
id2name = {}
#File.open('/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/2025-01-23-data-analysis-for-mnsc-and-plots/2023-08-18-results-stranded-IV-largeRNAFrac-RNAflow-COPY-relevant-data/deseq2_Avian_vs_Mock_filtered_padj_0.05_extended.csv').each do |line|
File.open('/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/2025-rnaseq-boxplots-for-paper/input-data/deseq2_Avian_vs_Mock_filtered_padj_0.05_extended.csv').each do |line|
    next if line.start_with?('ID,geneName')
    s = line.split(',')
    id2name[s[0]] = s[1]
end

Dir.glob("/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/2025-rnaseq-boxplots-for-paper/input-data/*_reps_tpms.tsv").each do |tsv|
    next if File.basename(tsv) == 'tpms.tsv' || File.basename(tsv) == 'pvals.tsv' 
    condition = File.basename(tsv, '_reps_tpms.tsv')
    conditions.push(condition)
    File.open(tsv, 'r').each do |line|
        next if line.start_with?('ID')
        s = line.split("\t")
        id = s[0].chomp
        gene_ids.push(id)
        tpm_rep1 = s[1].chomp; tpm_rep2 = s[2].chomp; tpm_rep3 = s[3].chomp
        #puts id+"\t"+tpm_rep1+"\t"+tpm_rep2+"\t"+tpm_rep3
        condition2tpm[condition][id] = [tpm_rep1, tpm_rep2, tpm_rep3]
    end
end

puts conditions
gene_ids.sort!.uniq!
puts gene_ids.size

output = File.open('/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/2025-rnaseq-boxplots-for-paper/input-data/tpms-human.tsv','w')
output << ['ID','Name','mock-rep1','mock-rep2','mock-rep3','swine-rep1','swine-rep2','swine-rep3','avian-rep1','avian-rep2','avian-rep3','reass-rep1','reass-rep2','reass-rep3'].join("\t") << "\n"

gene_ids.each do |gene_id|
    if id2name[gene_id]
        gene_name = id2name[gene_id]
    else
        gene_name = ''
    end
    output << gene_id << "\t" << gene_name
    ['mock', 'swine', 'avian', 'reass'].each do |condition|
        tpms = condition2tpm[condition][gene_id]
        output << "\t" << tpms.join("\t")
    end
    output << "\n"
end

output.close
