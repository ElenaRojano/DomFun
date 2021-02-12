#! /usr/bin/env ruby

##########################
# Rojano E. & Seoane P., Feb 2021
# Script to generate a report to consult the status of the 
# data in CAFA 3 and CATH for proteins function prediction
##########################

require 'optparse'

##########################
#METHODS
##########################

def load_hash(filename, domain_type='superfamilyID')
	uniprot_ACC_container = {}
	gene_ID_container = {}
	counter = 0
	File.open(filename).each do |line|
		line.chomp!
		if counter == 0
			counter += 1
			next
		end
		cath_data = line.split("\t", 9)
		uniprot_ACC = cath_data[0].gsub('"', '')
		gene_ID = cath_data[4].gsub('"', '')
		superfamily_ID = cath_data[5].gsub('"', '')
		funfam_ID = cath_data[6].gsub('"', '')
		if domain_type == 'superfamilyID'
			value = superfamily_ID
		else
			value = funfam_ID
		end	
		fill_hash(uniprot_ACC_container, uniprot_ACC, value)
		fill_hash(gene_ID_container, gene_ID, value)
	end
	return uniprot_ACC_container, gene_ID_container
end

def fill_hash(hash_name, key, value)
	query = hash_name[key]
	if query.nil?
		hash_name[key] = [value]
	else
		query << value unless query.include?(value) || value == '' #avoid empty domain records
	end
end

def load_array(filename)
	container = []
	File.open(filename).each do |line|
		line.chomp!
		value = line.split("\t")
		container << value
	end
	return container.uniq.flatten!
end

##########################
#OPT-PARSER
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"
  
  options[:cath_geneIDs] = nil
  opts.on("-a", "--cath_geneIDs PATH", "File from CATH with proteins and FF or SF") do |data|
    options[:cath_geneIDs] = data
  end

  options[:cafa3_targets] = nil
  opts.on("-b", "--cafa3_targets PATH", "CAFA 3 target proteins from all species") do |data|
    options[:cafa3_targets] = data
  end

  options[:cafa3_training] = nil
  opts.on("-c", "--cafa3_training PATH", "CAFA 3 training protein from all species") do |data|
    options[:cafa3_training] = data
  end

  options[:domain_category] = "superfamilyID"
  opts.on("-d", "--domain_category PATH", "Domain category. Please choose one: superfamilyID or funfamID" ) do |data|
    options[:domain_category] = data
  end

  options[:output_file] = 'output_stats.txt'
  opts.on("-o", "--output_file PATH", "Output file") do |data|
    options[:output_file] = data
  end

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!

##########################
#MAIN
##########################

stats = []

uniprot_ACC_container, gene_ID_container = load_hash(options[:cath_geneIDs], options[:domain_category])
#uniprot_ACC_container for CAFA3 training IDs
#gene_ID_container for CAFA3 testing IDs

cath_species = gene_ID_container.keys.map{|a| a.split('_').last}
cath_for_target_proteins = gene_ID_container.keys
cath_for_training_proteins = uniprot_ACC_container.keys

total_cath_geneIDs = cath_for_target_proteins.length
stats << ["Total genes in CATH:", total_cath_geneIDs]

cafa3_targets = load_array(options[:cafa3_targets])
total_cafa3_targets = cafa3_targets.length
stats << ["Total genes in CAFA3 (targets):", total_cafa3_targets]

matched_geneIDs = cath_for_target_proteins & cafa3_targets
total_target_genes_in_CATH = matched_geneIDs.uniq.length

stats << ['Genes with domains:', total_target_genes_in_CATH]
stats << ['Percentage of CAFA3 genes with at least one CATH domains:', total_target_genes_in_CATH.fdiv(total_cafa3_targets)*100]

unmatched_cafa3_targets = cafa3_targets - matched_geneIDs
total_cafa3_unmatched_targets = unmatched_cafa3_targets.length

cafa3_training = load_array(options[:cafa3_training])

total_cafa3_training = cafa3_training.length
stats << ["Total genes in CAFA3 (training):", total_cafa3_training]

matched_geneIDs_training = cath_for_training_proteins & cafa3_training
total_training_genes_in_CATH = matched_geneIDs_training.uniq.length
stats << ['Genes with domains:', total_training_genes_in_CATH]
stats << ['Percentage of CAFA3 genes with at least one CATH domains:', total_training_genes_in_CATH.fdiv(total_cafa3_training)*100]

File.open(options[:output_file], 'w') do |f|
	stats.each do |stat|
		f.puts stat.join("\t")
	end
end