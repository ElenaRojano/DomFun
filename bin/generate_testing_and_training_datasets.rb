#! /usr/bin/env ruby
##########################
# Rojano E. & Seoane P., July 2019
# Generate training and testing datasets from CAFA2 and UniProt data.
# This script translate identifiers and generate the file for the
# construction of tripartite network (domain-protein-FunSys)
##########################


REPORT_FOLDER=File.expand_path(File.join(File.dirname(__FILE__), '..', 'templates'))
ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'DomFun'))
require 'generalMethods.rb'
require 'optparse'
require 'csv'

##########################
#METHODS
##########################

def load_hash(filename, mode)
	container = {}
	File.open(filename).each do |line|
		line.chomp!
		key, value = line.split("\t") if mode == 'a'
		value, key = line.split("\t") if mode == 'b'
		container[key] = value
	end
	return container
end

def load_array(file)
	array = []
	File.open(file).each do |line|
		line.chomp!
		array << line
	end
	return array
end

def remove_testing_proteins_from_targets(testing, targets)
	targets.reject {|k,v| testing.include?(k)}
end

def load_protein_go_file(file)
	proteins_go_t0 = {}
	File.open(file).each do |line|
		line.chomp!
		proteinID, go_term = line.split("\t")
		query = proteins_go_t0[proteinID]
		if query.nil?
			proteins_go_t0[proteinID] = [go_term]
		else
			query << go_term unless proteins_go_t0[proteinID].include?(go_term)
		end
	end
	return proteins_go_t0
end

def translate_genename2proteinId(dictionary, genenames)
	training_proteinIDs = []
	untranslated_proteinIDs_list = []
	genenames.each do |a|
		unless dictionary[a].nil?
			training_proteinIDs << dictionary[a]
		else
			untranslated_proteinIDs_list << a
		end
	end
	return training_proteinIDs.uniq, untranslated_proteinIDs_list.uniq
end

def get_go_from_proteins(training_proteinIDs, proteins_go_t0)
	training_proteins = {}
	training_proteinIDs.each do |proteinID|
		go_terms = proteins_go_t0[proteinID]
		training_proteins[proteinID] = go_terms unless go_terms.nil?
	end
	return training_proteins
end

def build_tripartite_network(cath_data, training_proteins, filename)
	protein_go_pairs = []
	protein_domain_pairs = []
	training_proteins.each do |proteinId, go_terms|
		domains = cath_data[proteinId]
		unless domains.nil?
			domains.each do |domain|
				protein_domain_pairs << [domain, proteinId]
			end
			go_terms.each do |go_term|
				protein_go_pairs << [proteinId, go_term]
			end
		end
	end
	handler = File.open(filename, 'w')
	protein_go_pairs.each do |protein, go|
		handler.puts "#{go}\t#{protein}"
	end
	protein_domain_pairs.each do |domain, proteinId|
		handler.puts "#{domain}\t#{proteinId}"
	end
	handler.close
end

##########################
#OPT-PARSER
##########################
options = {}
OptionParser.new do |opts|
	opts.banner = "Usage: #{__FILE__} [options]"

	options[:uniprotACC_file] = nil
		opts.on("-a", "--uniprotACC_file PATH", "Uniprot ACC identifiers file") do |data|
		options[:uniprotACC_file] = data
	end

	options[:cafa2_identifiers] = nil
		opts.on("-c", "--cafa2_identifiers PATH", "CAFA2 identifiers file, including UniProt ACC") do |data|
		options[:cafa2_identifiers] = data
	end

	options[:cath_file] = nil
		opts.on("-C", "--cath_file PATH", "CATH file to get protein domains") do |data|
		options[:cath_file] = data
	end

	
	options[:targets_file] = nil
		opts.on("-g", "--targets_file PATH", "CAFA2 and genename ids from FASTA files") do |data|
		options[:targets_file] = data
	end

	options[:output_network] = 'tripartite_network.txt'
		opts.on("-o", "--output_network PATH", "Output tripartite network from CAFA information") do |data|
		options[:output_network] = data
	end

	options[:proteins_go] = nil
		opts.on("-p", "--proteins_go PATH", "Protein identifiers with GO terms at time 0") do |data|
		options[:proteins_go] = data
	end

	options[:testing_file] = nil
		opts.on("-t", "--testing_file PATH", "Testing file with CAFA2 accession identifiers") do |data|
		options[:testing_file] = data
	end

	opts.on_tail("-h", "--help", "Tool information") do
		puts opts
		exit
	end

end.parse!

##########################
#MAIN
##########################

testing_proteins_CAFA2_ids = load_array(options[:testing_file])
targets_for_training = load_hash(options[:targets_file], 'a')

# Eliminate testing proteins from target (training) to avoid redundancy:
remove_testing_proteins_from_targets(testing_proteins_CAFA2_ids, targets_for_training)

gene_term_dictionary = load_hash(options[:uniprotACC_file], 'b')

training_proteinIDs, untranslated_proteinIDs_list = translate_genename2proteinId(gene_term_dictionary, targets_for_training.values)

proteins_go_t0 = load_protein_go_file(options[:proteins_go])
training_proteins = get_go_from_proteins(training_proteinIDs, proteins_go_t0)

cath_data, protein2gene, gene2proteins, cath_proteins_number = load_cath_data(options[:cath_file], 'funfamID')

build_tripartite_network(cath_data, training_proteins, options[:output_network])