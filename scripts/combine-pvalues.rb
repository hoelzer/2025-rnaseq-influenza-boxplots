id2adjp = {}
gene_ids = []
Dir.glob('/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/2025-01-23-data-analysis-for-mnsc-and-plots/2023-08-18-results-stranded-IV-largeRNAFrac-RNAflow-COPY-relevant-data/*full*.csv').each do |csv|
    basename = File.basename(csv,'_full_extended.csv').split('deseq2_')[1]
    condition1 = basename.split('_vs_')[0]
    condition2 = basename.split('_vs_')[1].chomp
    if condition2 == 'Mock'
        condition2 = condition1 # Mock becomes Avian
        condition1 = 'Mock'
    end
    condition1 = 'mock' if condition1 == 'Mock'
    condition2 = 'avian' if condition2 == 'Avian'
    condition1 = 'avian' if condition1 == 'Avian'
    condition2 = 'swine' if condition2 == 'Swine'
    condition1 = 'swine' if condition1 == 'Swine'
    condition2 = 'reass' if condition2 == 'Reassortant'
    condition1 = 'reass' if condition1 == 'Reassortant'

    puts "#{condition1}-#{condition2}"
    id2adjp["#{condition1}-#{condition2}"] = {}

    File.open(csv,'r').each do |line|
        next if line.start_with?('ID,geneName')
        s = line.split(',')
        gene_id = s[0]
        adjp = s[11]
        gene_ids.push(gene_id)
        id2adjp["#{condition1}-#{condition2}"][gene_id] = adjp
    end
end

gene_ids.sort!.uniq!
puts "#{gene_ids.size} Gene IDs found."

#output = File.open('pvals.tsv','w')
output = File.open('pvals-all-human-comparisons.tsv','w')
output << "ID"
id2adjp.keys.each do |comparison|
    output << "\t#{comparison}"
end
output << "\n"

gene_ids.each do |gene_id|
    output << "#{gene_id}"
    id2adjp.each do |comparison, id2adjp_value|
        if id2adjp_value[gene_id]
            output << "\t#{id2adjp_value[gene_id]}"
        else
            output << "\t0.05" # fill up missing pvalues. Should not happen bc using now the full tables from DESeq2
        end
    end
    output << "\n"
end
output.close
