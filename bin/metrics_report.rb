#! /usr/bin/env ruby

##########################
# Rojano E. & Seoane P., Feb 2021
# Script to generate a report to consult the status of the 
# data in CAFA and CATH for proteins function prediction
##########################

ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'DomFun'))
require 'generalMethods'
require 'optparse'
require 'csv'

##########################
#METHODS
##########################

def load_data(input_file)
	metadata = []
	File.open(input_file).each do |line|
		line.chomp!
		metadata << line.split("\t")
	end
	return metadata
end

def load_tab_file(training_proteins)
	tab_data = []
	File.open(training_proteins).each do |line|
		line.chomp!
		tab_data << line.split("\t")
	end
	return tab_data
end

def calculate_training_proteins_stats(training_proteins, stats, prefix)
	stats[prefix + 'number'] = training_proteins.map{|a| a.first}.uniq.length
	stats[prefix + 'all_go_annots'] = training_proteins.map{|a| a[1]}.uniq.length
	go_subontologies = {'C' => [], 'F' => [], 'P' => []}
	training_proteins.each do |protID, goID, goType|
		go_subontologies[goType] << goID
	end
	stats[prefix + 'total_gocc'] = go_subontologies['C'].uniq.length
	stats[prefix + 'total_gomf'] = go_subontologies['F'].uniq.length
	stats[prefix + 'total_gobp'] = go_subontologies['P'].uniq.length
	stats[prefix + 'total_protein-gocc_relations'] = go_subontologies['C'].length
	stats[prefix + 'total_protein-gocc_gomf'] = go_subontologies['F'].length
	stats[prefix + 'total_protein-gocc_gobp'] = go_subontologies['P'].length
end

def calculate_testing_proteins_stats(testing_proteins, stats, prefix)
	stats[prefix + 'number'] = testing_proteins.map{|a| a.first}.uniq.length
	proteins_by_organisms = testing_proteins.map{|a| a.first.split('_').last}
	list_of_proteins_by_organisms = proteins_by_organisms.group_by{|e| e}.map{|k, v| [k, v.length]}.to_h
	list_of_proteins_by_organisms.each do |organism, count|
		stats[prefix + organism + '_proteins'] = count if organism == 'HUMAN'
	end
end


def calculate_cath_stats(cath_info_domains, stats, prefix)
	stats[prefix + 'number'] = cath_info_domains.values.uniq.length
	stats[prefix + 'protein_relations'] = cath_info_domains.values.map{|v| v.uniq.length}.inject { |sum, n| sum + n } 
end


def get_input_cafa_data(path, stats, id)		
	training_proteins = File.join(path, 'training_proteins.txt')
	training_proteins = load_tab_file(training_proteins)
	prefix = id.upcase + '==training_proteins--'
	calculate_training_proteins_stats(training_proteins, stats, prefix)

	human_training_proteins = File.join(path, 'training_proteins_human.txt')
	training_proteins = load_tab_file(human_training_proteins)
	prefix = id.upcase + '==HUMAN_training_proteins--'
	calculate_training_proteins_stats(training_proteins, stats, prefix)

	testing_proteins = File.join(path, 'testing_proteins.txt')
	testing_proteins = load_tab_file(testing_proteins)
	prefix = id.upcase + '==testing_proteins--'
	calculate_testing_proteins_stats(testing_proteins, stats, prefix)
end

def get_input_CATH_data(path, stats, id)
	#cath_info_SF == {:ProtID => ["SFid1", "SFid2"]}
	#protein2gene_dict_SF == {:ProtID => "GeneName"}
	#cath_proteins_number == 38581 (TOTAL PROTEINS WITH DOMAINS ANNOTATIONS)
		
	cath_info_SF, protein2gene_dict_SF, cath_proteins_number = load_cath_data(path, 'superfamilyID', whitelist=nil, dictionary_key='gene_name')
	prefix = id.upcase + '==superfamily_domains--'
	calculate_cath_stats(cath_info_SF, stats, prefix)	
	
	cath_info_FF, protein2gene_dict_FF, cath_proteins_number = load_cath_data(path, 'funfamID', whitelist=nil, dictionary_key='gene_name')
	prefix = id.upcase + '==funfam_domains--'
	calculate_cath_stats(cath_info_FF, stats, prefix)
	
	stats[id.upcase + '==proteins_number'] = cath_proteins_number #only once
end

def get_input_network_data(path, stats, id)
	Dir.glob(path).each do |input_file|
		if input_file.include?('superfamilyID')
			if input_file.include?('GOMF')
				prefix = id.upcase + '==superfamilyGOMFNetwork--'
				calculate_network_stats(input_file, stats, prefix, id)
			elsif input_file.include?('GOBP')
				prefix = id.upcase + '==superfamilyGOBPNetwork--'
				calculate_network_stats(input_file, stats, prefix, id)
			elsif input_file.include?('GOCC')
				prefix = id.upcase + '==superfamilyGOCCNetwork--'
				calculate_network_stats(input_file, stats, prefix, id)
			end
		end	
	end
end

def calculate_network_stats(input_file, stats, prefix, id)
	network_info = {}
	File.open(input_file).each do |line|
		line.chomp!
		identifier, proteinID = line.split("\t")
		query = network_info[proteinID]
		if query.nil?
			network_info[proteinID] = [identifier]
		else
			query << identifier
		end
	end
	stats[id.upcase + '==uniq_proteins_number'] = network_info.keys.flatten.uniq.length
	stats[id.upcase + '==GO-proteins_relations'] = network_info.values.flatten.select {|el| el.include?('GO') }.length
	stats[id.upcase + '==domain-proteins_relations'] = network_info.values.flatten.select {|el| el.include?('.') }.length
	stats[id.upcase + '==total_domain/GO-protein_relations'] = network_info.values.flatten.length
end

def get_input_associations_data(path, stats, id)

end


##########################
#OPT-PARSER
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"
  
  options[:input_file] = nil
  opts.on("-i", "--input_file PATH", "Input file with tags and paths to analyze") do |data|
    options[:input_file] = data
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


tag_path_info = load_data(options[:input_file])
stats = {}
tag_path_info.each do |id, path|
	if id.include?('input_cafa')
		get_input_cafa_data(path, stats, id)
	elsif id.include?('input_cath')
		get_input_CATH_data(path, stats, id)
	elsif id.include?('input_networks')
		get_input_network_data(path, stats, id)
	elsif id.include?('input_associations')
		get_input_associations_data(path, stats, id)
	end
end

File.open(options[:output_file], 'w') do |f|
	stats.each do |stat, values|
		f.puts "#{stat}\t#{values}"
	end
end

Process.exit