############################################## DESeq2 normalized counts, old code

#output = File.open('/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/TPMs/counts-segments.tsv','w')
#output << ['ID','Name','avian-rep1','avian-rep2','avian-rep3','mock-rep1','mock-rep2','mock-rep3','reass-rep1','reass-rep2','reass-rep3','swine-rep1','swine-rep2','swine-rep3'].join("\t") << "\n"
#
#File.open('/Users/martin/projects/2024-04-16-manual-DESeq-Agustina/2025-04-08/results/data/counts/normalized_counts.csv','r').each do |line|
#    next if line.start_with?('"","avian')
#    s = line.split(",")
#    id = s[0].gsub('"','').chomp
#    name = id.sub('gene-','')
#    next if id.start_with?('ENSG')
#    counts = s[1..s.length]
#    output << id << "\t" << name << "\t" << counts.join("\t")
#end
#
#output.close


######################################################## TPM counts 

output = File.open('/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/2025-rnaseq-boxplots-for-pape/input-data/counts-tpm-segments.tsv','w')
output << ['ID','Name','avian-rep1','avian-rep2','avian-rep3','mock-rep1','mock-rep2','mock-rep3','reass-rep1','reass-rep2','reass-rep3','swine-rep1','swine-rep2','swine-rep3'].join("\t") << "\n"

condition2tsv = {}
gene_ids = []

Dir.glob("/Users/martin/projects/2025-03-13-Influenza-RNASeq-Agustina/2025-rnaseq-boxplots-for-paper/input-data/segments-strand2/*.tsv").each do |tsv|
    condition = File.basename(tsv, '.counts.tpm.tsv').gsub('_','-').gsub('reassortant','reass')
    puts condition
    condition2tsv[condition] = {}
    File.open(tsv,'r').each do |line|
        next if line.start_with?('Geneid')
        s = line.split("\t")
        id = s[0].chomp
        name = id.sub('gene-','')
        next if id.start_with?('ENSG')
        count = s[7].chomp
        condition2tsv[condition][name] = count
        gene_ids.push(name)
    end
end

gene_ids.sort!.uniq!

gene_ids.each do |gene_name|
    output << "gene-" << gene_name << "\t" << gene_name
    ['avian-rep1','avian-rep2','avian-rep3','mock-rep1','mock-rep2','mock-rep3','reass-rep1','reass-rep2','reass-rep3','swine-rep1','swine-rep2','swine-rep3'].each do |condition|
        output << "\t" << condition2tsv[condition][gene_name] 
    end
    output << "\n"
end

output.close
