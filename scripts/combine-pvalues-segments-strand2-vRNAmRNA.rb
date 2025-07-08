id2adjp = {}
gene_ids = []
# cp /Users/martin/projects/2024-04-16-manual-DESeq-Agustina/2025-05-19-deseq2-for-vRNAmRNA/results/*vs*/results/*full*.csv ~/projects/2025-03-13-Influenza-RNASeq-Agustina/2025-rnaseq-boxplots-for-paper/input-data/segments-strand2-vRNAmRNA
Dir.glob('/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/2025-rnaseq-boxplots-for-paper/input-data/segments-strand2-vRNAmRNA/*full*.csv').each do |csv|
    basename = File.basename(csv,'_full.csv').split('deseq2_')[1]
    condition1 = basename.split('_vs_')[0]
    condition2 = basename.split('_vs_')[1].chomp

    condition1 = 'mock' if condition1 == 'Mock'
    condition2 = 'avian' if condition2 == 'Avian'
    condition2 = 'swine' if condition2 == 'Swine'
    condition2 = 'reass' if condition2 == 'Reassortant'
    condition1 = 'avian' if condition2 == 'Avian'
    condition1 = 'swine' if condition2 == 'Swine'
    condition1 = 'reass' if condition2 == 'Reassortant'

    puts "#{condition1}-#{condition2}"
    id2adjp["#{condition1}-#{condition2}"] = {}

    File.open(csv,'r').each do |line|
        next if line.start_with?('"","baseMean')
        s = line.split(',')
        gene_id = s[0].gsub('"','')
        next if gene_id.start_with?('ENSG')
        adjp = s[5].chomp
        gene_ids.push(gene_id)
        id2adjp["#{condition1}-#{condition2}"][gene_id] = adjp
    end
end

gene_ids.sort!.uniq!
puts "#{gene_ids.size} Gene IDs found from the Influenza segments."

output = File.open('input-data/pvals-segments-strand2-vRNAmRNA.tsv','w')
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
            puts "ATTENTON: dummy pvalue bc missing"
            output << "\t0.05" # fill up missing pvalues. Should not happen bc using now the full tables from DESeq2
        end
    end
    output << "\n"
end
output.close
