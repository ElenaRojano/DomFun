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

def open_file_save_array(input_file, output_array)
	File.open(input_file).each do |line|
		line.chomp!
		output_array << line.split("\t")
	end
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
				prefix = id.upcase + '==superfamily_GOMF_network--'
				calculate_network_stats(input_file, stats, prefix, id)
			elsif input_file.include?('GOBP')
				prefix = id.upcase + '==superfamily_GOBP_network--'
				calculate_network_stats(input_file, stats, prefix, id)
			elsif input_file.include?('GOCC')
				prefix = id.upcase + '==superfamily_GOCC_network--'
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
	stats[prefix + 'uniq_proteins_number'] = network_info.keys.flatten.uniq.length
	stats[prefix + 'uniq_domains_number'] = network_info.values.flatten.select{|a| a.include?('.')}.uniq.length
	stats[prefix + 'uniq_GO_annotations_number'] = network_info.values.flatten.select{|a| a.include?('GO')}.uniq.length
	stats[prefix + 'GO-proteins_relations'] = network_info.values.flatten.select {|el| el.include?('GO') }.length
	stats[prefix + 'domain-proteins_relations'] = network_info.values.flatten.select {|el| el.include?('.') }.length
	stats[prefix + 'total_domain/GO-protein_relations'] = network_info.values.flatten.length
end

#Note: associations, predictions and normalized predictions
#load file are under the same method.
#Stats calculation has been joined in the same method aw.
def get_input_assoc_predictions_data(path, stats, id, info_string)
	assoc_methods = ['jaccard', 'simpson', 'hypergeometric', 'pcc']
	if info_string == '_assocs--'
		mode = 'associations'
	elsif info_string == '_preds--'
		mode = 'predictions'
	elsif info_string == '_norm_preds--'
		mode = 'norm_predictions'
	end		
	Dir.glob(path).each do |input_file|
		assoc_methods.each do |assoc_method|		
			if input_file.include?('superfamilyID')
				if input_file.include?('GOMF') && input_file.include?(assoc_method)
					prefix = id.upcase + '==superfamily_GOMF_' + assoc_method + info_string
					calculate_stats(input_file, stats, prefix, id, mode)
				elsif input_file.include?('GOBP') && input_file.include?(assoc_method)
					prefix = id.upcase + '==superfamily_GOBP_' + assoc_method + info_string
 					calculate_stats(input_file, stats, prefix, id, mode)
				elsif input_file.include?('GOCC') && input_file.include?(assoc_method)
					prefix = id.upcase + '==superfamily_GOCC' + assoc_method + info_string
					calculate_stats(input_file, stats, prefix, id, mode)
				end
			end
		end	
	end
end

def calculate_stats(input_file, stats, prefix, id, mode)
	file_info = []
	open_file_save_array(input_file, file_info)
	if mode == 'associations'
		stats[prefix + 'GO-domain_relations_number'] = file_info.length
		stats[prefix + 'GO_number'] = file_info.map{|a| a.first}.uniq.length
		stats[prefix + 'domain_number'] = file_info.map{|a| a[1]}.uniq.length
		stats[prefix + 'max_assoc_val'] = file_info.map{|a| a.last}.max
		stats[prefix + 'min_assoc_val'] = file_info.map{|a| a.last}.min
		assoc_vals = file_info.map{|a| a.last.to_f}
		stats[prefix + 'ave_assoc_val'] = assoc_vals.inject { |sum, n| sum + n }.fdiv(assoc_vals.length)
	elsif mode == 'predictions'
		stats[prefix + 'predicted_proteins'] = file_info.map{|a| a.first}.uniq.length
		stats[prefix + 'total_GOs_predicted'] = file_info.map{|a| a[2]}.uniq.length
		stats[prefix + 'total_domains'] = file_info.map{|a| a[1].split(',')}.flatten.uniq.length
		stats[prefix + 'max_pred_val'] = file_info.map{|a| a.last}.max
		stats[prefix + 'min_pred_val'] = file_info.map{|a| a.last}.min
		pred_vals = file_info.map{|a| a.last.to_f}
		stats[prefix + 'ave_pred_val'] = pred_vals.inject { |sum, n| sum + n }.fdiv(pred_vals.length)
	elsif mode == 'norm_predictions'
		#Norm predictions max value are wrong (gets scientific annotation like: 9.101087297347377e-06 as max)
		#See how to correct it
		stats[prefix + 'predicted_proteins_after_norm'] = file_info.map{|a| a.first}.uniq.length
		stats[prefix + 'total_GOs_predicted_after_norm'] = file_info.map{|a| a[2]}.uniq.length
		stats[prefix + 'total_domains_after_norm'] = file_info.map{|a| a[1].split(',')}.flatten.uniq.length
		stats[prefix + 'max_norm_pred_val'] = file_info.map{|a| a.last}.max
		stats[prefix + 'min_norm_pred_val'] = file_info.map{|a| a.last}.min
		pred_vals = file_info.map{|a| a.last.to_f}
		stats[prefix + 'ave_norm_pred_val'] = pred_vals.inject { |sum, n| sum + n }.fdiv(pred_vals.length)
	end
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
loaded_data ={}
tag_path_info.each do |id, path|
	if id.include?('input_cafa')
		get_input_cafa_data(path, stats, id)
	elsif id.include?('input_cath')
		get_input_CATH_data(path, stats, id)
	elsif id.include?('input_networks')
		get_input_network_data(path, stats, id)
	elsif id.include?('input_associations')
		info_string = '_assocs--'
		get_input_assoc_predictions_data(path, stats, id, info_string)
	elsif id.include?('input_predictions')
		info_string = '_preds--'
		get_input_assoc_predictions_data(path, stats, id, info_string)
	elsif id.include?('input_normpreds')
		info_string = info_string = '_norm_preds--'
		get_input_assoc_predictions_data(path, stats, id, info_string)
	#TODO: calculate how many proteins have CATH domains
	end
end



File.open(options[:output_file], 'w') do |f|
	stats.each do |stat, values|
		f.puts "#{stat}\t#{values}"
	end
end

Process.exit