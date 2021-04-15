#! /usr/bin/env ruby

##########################
# Rojano E. & Seoane P., Feb 2021
# Script to generate a report to consult the status of the 
# data in CAFA and CATH for proteins function prediction
##########################

REPORT_FOLDER=File.expand_path(File.join(File.dirname(__FILE__), '..', 'templates'))
ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'DomFun'))
require 'generalMethods'
require 'colorize'
require 'optparse'
require 'csv'
require 'report_html'
require 'json'

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

def load_tab_file(file)
	tab_data = []
	File.open(file).each do |line|
		line.chomp!
		tab_data << line.split("\t")
	end
	return tab_data
end

def calculate_training_proteins_stats(training_proteins, testing_training_stats)
	testing_training_stats['training'] = {'general_stats' => {}, 'GOMF' => {}, 'GOBP' => {}, 'GOCC' => {}}
	go_subontologies = {'C' => [], 'F' => [], 'P' => []}
	all_go_terms = []
	annotated_proteins = []
	training_proteins.each do |protID, goID, goType|
		go_subontologies[goType] << goID
		all_go_terms << goID
		annotated_proteins << [protID, goType]
	end

	testing_training_stats['training']['general_stats']['terms'] = all_go_terms.uniq.length
	testing_training_stats['training']['general_stats']['proteins'] = annotated_proteins.map{|a| a.first}.uniq.length

	gomf_annotated_proteins = annotated_proteins.select{|a| a.last == 'F'}
	gobp_annotated_proteins = annotated_proteins.select{|a| a.last == 'P'}
	gocc_annotated_proteins = annotated_proteins.select{|a| a.last == 'C'}

	testing_training_stats['training']['GOMF']['proteins'] = gomf_annotated_proteins.uniq.length
	testing_training_stats['training']['GOBP']['proteins'] = gobp_annotated_proteins.uniq.length
	testing_training_stats['training']['GOCC']['proteins'] = gocc_annotated_proteins.uniq.length
	

	testing_training_stats['training']['GOMF']['terms'] = go_subontologies['F'].uniq.length
	testing_training_stats['training']['GOBP']['terms'] = go_subontologies['P'].uniq.length
	testing_training_stats['training']['GOCC']['terms'] = go_subontologies['C'].uniq.length

	#TODO: protein_GOMF_relations and proteins_annotated_in give the same value.
	#Select which one to add.
	testing_training_stats['training']['GOMF']['protein-term'] = go_subontologies['F'].length
	testing_training_stats['training']['GOBP']['protein-term'] = go_subontologies['P'].length
	testing_training_stats['training']['GOCC']['protein-term'] = go_subontologies['C'].length

	return gomf_annotated_proteins, gobp_annotated_proteins, gocc_annotated_proteins
end

def calculate_testing_proteins_stats(testing_proteins, testing_training_stats, human_testing_training_stats, geneAccession_protID_dictionary, testing_storage)

	testing_training_stats['testing'] = {}
	human_testing_training_stats['testing'] = {}
	human_proteins = testing_proteins.select{|prot| prot.include?('_HUMAN')}
	
	all_testing_proteins_translated, all_testing_proteins_untranslated = translate_protID2geneName(testing_proteins, geneAccession_protID_dictionary)
	human_testing_proteins_translated, human_testing_proteins_untranslated = translate_protID2geneName(human_proteins, geneAccession_protID_dictionary)

	proteins_by_organisms = testing_proteins.map{|a| a.split('_').last}
	list_of_proteins_by_organisms = proteins_by_organisms.group_by{|e| e}.map{|k, v| [k, v.length]}.to_h
	
	human_testing_training_stats['testing']['proteins_geneID'] = list_of_proteins_by_organisms['HUMAN']
	human_testing_training_stats['testing']['proteins_UniProtID_translated'] = human_testing_proteins_translated.length
	human_testing_training_stats['testing']['proteins_UniProtID_untranslated'] = human_testing_proteins_untranslated.length
	
	testing_training_stats['testing']['proteins_geneID'] = list_of_proteins_by_organisms.values.inject{|a, b| a + b}
	testing_training_stats['testing']['proteins_UniProtID_translated'] = all_testing_proteins_translated.length
	testing_training_stats['testing']['proteins_UniProtID_untranslated'] = all_testing_proteins_untranslated.length

	testing_storage['human'] = human_testing_proteins_translated
	testing_storage['all'] = all_testing_proteins_translated

end


def calculate_cath_stats(cath_info_domains, id, domain_type, domains_stats)
	#Leave this method in case you need to add more stats
	domains_stats[domain_type]['number'] = cath_info_domains.values.uniq.length
end

def load_geneAccession_protID_dictionary(dictionary_file, geneAccession_protID_dictionary)
	File.open(dictionary_file).each do |line|
		line.chomp!
		protID, geneAccessionID = line.split("\t")
		geneAccession_protID_dictionary[geneAccessionID] = protID
	end
end

def translate_protID2geneName(proteinIDs, geneAccession_protID_dictionary)
	testing_proteins_translated = []
	testing_proteins_untranslated = []
	proteinIDs.each do |geneAccessionID|
		proteinID = geneAccession_protID_dictionary[geneAccessionID]
		unless proteinID.nil?
			testing_proteins_translated << proteinID
		else
			testing_proteins_untranslated << geneAccessionID
		end
	end	
	return testing_proteins_translated, testing_proteins_untranslated
end

def load_network_files(input_file)
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
	return network_info
end

def get_input_cafa_data(path, id, geneAccession_protID_dictionary, stats_complex, training_storage, testing_storage)
	all_testing_training_stats = {'testing' => {}, 'training' => {}}
	human_testing_training_stats = {'testing' => {}, 'training' => {}}
	training_storage = {'human' => {'GOMF' => [], 'GOBP' => [], 'GOCC' => []}, 'all' => {'GOMF' => [], 'GOBP' => [], 'GOCC' => []}}
	testing_storage = {'human' => [], 'all' => []}

	training_proteins = File.join(path, 'training_proteins.txt')
	training_proteins = load_tab_file(training_proteins)
	
	all_gomf_annotated_proteins, all_gobp_annotated_proteins, all_gocc_annotated_proteins = calculate_training_proteins_stats(training_proteins, all_testing_training_stats)
	
	#Info for combined results:
	training_storage['all']['GOMF'] = all_gomf_annotated_proteins
	training_storage['all']['GOBP'] = all_gobp_annotated_proteins
	training_storage['all']['GOCC'] = all_gocc_annotated_proteins
	
	human_training_proteins = File.join(path, 'training_proteins_human.txt')
	human_training_proteins = load_tab_file(human_training_proteins)
	human_gomf_annotated_proteins, human_gobp_annotated_proteins, human_gocc_annotated_proteins = calculate_training_proteins_stats(human_training_proteins, human_testing_training_stats)

	training_storage['human']['GOMF'] = human_gomf_annotated_proteins
	training_storage['human']['GOBP'] = human_gobp_annotated_proteins
	training_storage['human']['GOCC'] = human_gocc_annotated_proteins

	testing_proteins = File.join(path, 'testing_proteins.txt')
	testing_proteins = load_tab_file(testing_proteins)
	testing_proteins.flatten!
	calculate_testing_proteins_stats(testing_proteins, all_testing_training_stats, human_testing_training_stats, geneAccession_protID_dictionary, testing_storage)
	stats_complex['cafa'] = {'human' => human_testing_training_stats, 'all' => all_testing_training_stats}

end

def get_input_CATH_data(path, id, stats_complex)
	#cath_info_SF == {:ProtID => ["SFid1", "SFid2"]}
	#protein2gene_dict_SF == {:ProtID => "GeneName"}
	#cath_proteins_number == 38581 (TOTAL PROTEINS WITH DOMAINS ANNOTATIONS)
	
	superfamily_ID = 'superfamilyID'
	funfam_ID = 'funfamID'
	domains_stats = {superfamily_ID => {}, funfam_ID => {}}

	cath_info_SF, protein2gene_dict, cath_proteins_number = load_cath_data(path, superfamily_ID, whitelist=nil, dictionary_key='gene_name')	
	calculate_cath_stats(cath_info_SF, id, superfamily_ID, domains_stats)	

	cath_info_FF, protein2gene_dict, cath_proteins_number = load_cath_data(path, funfam_ID, whitelist=nil, dictionary_key='gene_name')
	calculate_cath_stats(cath_info_FF, id, funfam_ID, domains_stats)
	
	domains_stats['CATH_proteins_number'] = cath_proteins_number

	stats_complex['domains'] = domains_stats

end

def get_input_network_data(path, id, stats_complex)
	#network_data = {'human' => {}, 'all' => {}}
	if id.include?('human')
		orgtag = 'human'
	elsif id.include?('all')
		orgtag = 'all'
	end
	network_data_complex = {orgtag => {}}
	Dir.glob(path).each do |input_file|
		#STDERR.puts input_file
		if input_file.include?('superfamilyID')
			domtag = 'superfamilyID'
		elsif input_file.include?('funfamID')
			domtag = 'funfamID'
		end
		if input_file.include?('GOMF')
			annot = 'GOMF'
		elsif input_file.include?('GOBP')
			annot = 'GOBP'
		elsif input_file.include?('GOCC')
			annot = 'GOCC'
		end
		network_info = load_network_files(input_file)
		calculate_network_stats(network_data_complex[orgtag], domtag, orgtag, annot, network_info)
	end
	stats_complex['network'] = network_data_complex
end

def calculate_network_stats(network_data_complex, domtag, orgtag, annot, network_info)
	query_domtag = attrib_to_hash(network_data_complex, domtag, {annot => {}})
	query_annot = attrib_to_hash(query_domtag, annot)
	values_flatten = network_info.values.flatten
	query_annot['proteins'] = network_info.keys.flatten.uniq.length
	query_annot['domains'] = values_flatten.select{|a| a.include?('.')}.uniq.length
	query_annot['terms'] = values_flatten.select{|a| a.include?('GO')}.uniq.length
	query_annot['term-protein'] = values_flatten.select {|el| el.include?('GO') }.length
	query_annot['domain-protein'] = values_flatten.select {|el| el.include?('.') }.length
end

def attrib_to_hash(hash, key, value={})
	query = hash[key]
	if query.nil?
		query = value
		hash[key] = query
	end
	return query
end

#Note: associations, predictions and normalized predictions
#load file are under the same method.
#Stats calculation has been joined in the same method aw.
def get_input_assoc_predictions_data(path, id, mode, stats_complex)
	assoc_methods = ['jaccard', 'simpson', 'hypergeometric', 'pcc']
	if id.include?('human')
		orgtag = 'human'
	elsif id.include?('all')
		orgtag = 'all'
	end
	
	Dir.glob(path).each do |input_file|
		stats = {}
		assoc_method = assoc_methods.select{|am| input_file.include?(am)}.first	
		if input_file.include?('superfamilyID')
			domtag = 'superfamilyID'
		elsif input_file.include?('funfamID')
			domtag = 'funfamID'
		end
		if input_file.include?('GOMF')
			annot = 'GOMF'
		elsif input_file.include?('GOBP')
			annot = 'GOBP'
		elsif input_file.include?('GOCC')
			annot = 'GOCC'
		end
		file_info = load_tab_file(input_file)
		calculate_stats(file_info, mode, domtag, orgtag, annot, assoc_method, stats)
		stats_complex[mode] = stats
	end
end

#def calculate_stats(input_file, stats, prefix, id, mode, loaded_data)
def calculate_stats(file_info, mode, domtag, orgtag, annot, assoc_method, annotation_info_complex)
	query_orgtag = attrib_to_hash(annotation_info_complex, orgtag, {domtag => {}})
	query_domtag = attrib_to_hash(query_orgtag, domtag, {annot => {}})
	query_annot = attrib_to_hash(query_domtag, annot, {assoc_method => {}})
	query_assoc = attrib_to_hash(query_annot, assoc_method)

	if mode == 'assocs'
		query_assoc['term-domain'] = file_info.length
		query_assoc['terms'] = file_info.map{|a| a.first}.uniq.length
		query_assoc['domains'] = file_info.map{|a| a[1]}.uniq.length
		query_assoc['max_assoc_val'] = file_info.map{|a| a.last}.max
		query_assoc['min_assoc_val'] = file_info.map{|a| a.last}.min
		assoc_vals = file_info.map{|a| a.last.to_f}
		query_assoc['ave_assoc_val'] = assoc_vals.inject { |sum, n| sum + n }.fdiv(assoc_vals.length)
	elsif mode == 'preds' || mode == 'norm_preds'
		query_assoc['proteins'] = file_info.map{|a| a.first}.uniq.length
		query_assoc['terms'] = file_info.map{|a| a[2]}.uniq.length
		query_assoc['domains'] = file_info.map{|a| a[1].split(',')}.flatten.uniq.length
		query_assoc['max_pred_val'] = file_info.map{|a| a.last}.max
		query_assoc['min_pred_val'] = file_info.map{|a| a.last}.min
		pred_vals = file_info.map{|a| a.last.to_f}
		query_assoc['ave_pred_val'] = pred_vals.inject { |sum, n| sum + n }.fdiv(pred_vals.length)
	end
end

#TODO: CHANGE THIS METHOD TO THE NEW SYSTEM:
def get_combined_stats(loaded_data, stats, path, geneAccession_protID_dictionary)
	cath_tag = nil
	if path.include?('old')
		cath_tag = 'OLD'
	elsif path.include?('cur')
		cath_tag = 'CUR'
	end

	if path.include?('human')
		# 14/04 -> continuar para solucionar el problema de las relaciones go-prot del training
		# STDERR.puts loaded_data['human_training_proteins'].inspect
		# Process.exit
		training_proteins = loaded_data['human_training_proteins'].map{|a| a.first}.uniq
		testing_proteins = loaded_data['human_testing_proteins_translated'].uniq
		testing_proteins_untranslated = loaded_data['human_testing_proteins_untranslated'].uniq
		calculate_combined_stats(stats, path, loaded_data, 'HUMAN', training_proteins, testing_proteins, testing_proteins_untranslated, cath_tag)
	elsif path.include?('all')
		training_proteins = loaded_data['training_proteins'].map{|a| a.first}.uniq
		testing_proteins = loaded_data['testing_proteins_translated'].uniq
		testing_proteins_untranslated = loaded_data['testing_proteins_untranslated'].uniq
		calculate_combined_stats(stats, path, loaded_data, 'ALL', training_proteins, testing_proteins, testing_proteins_untranslated, cath_tag)
	end
end

def calculate_combined_stats(stats, path, loaded_data, organism, training_proteins, testing_proteins, testing_proteins_untranslated, cath_tag)
go_subontologies = ['GOMF', 'GOBP', 'GOCC']
	go_subontologies.each do |go_type|
		if path.include?('superfamilyID') && path.include?(go_type)
			cath_proteins = loaded_data['cath_info_SF'].keys.map{|a| a.to_s}
			stats["COMBINED_RESULTS" + '_' + cath_tag + '_' + organism + "==superfamilyID_" + go_type + '--intersection_training_cath_SF'] = (training_proteins & cath_proteins).length
			stats["COMBINED_RESULTS" + '_' + cath_tag + '_' + organism + "==superfamilyID_" + go_type + '--intersection_training_cath_SF (%)'] = (training_proteins & cath_proteins).length.fdiv(training_proteins.length)*100
			stats["COMBINED_RESULTS" + '_' + cath_tag + '_' + organism + "==superfamilyID_" + go_type + '--intersection_testing_cath_SF'] = (testing_proteins & cath_proteins).length
			stats["COMBINED_RESULTS" + '_' + cath_tag + '_' + organism + "==superfamilyID_" + go_type + '--intersection_testing_cath_SF (%)'] = (testing_proteins & cath_proteins).length.fdiv(testing_proteins.length)*100
		elsif path.include?('funfamID') && path.include?(go_type)
			cath_proteins = loaded_data['cath_info_FF'].keys.map{|a| a.to_s}
			stats["COMBINED_RESULTS" + '_' + cath_tag + '_' + organism + "==funfamID_" + go_type +'--intersection_training_cath_FF'] = (training_proteins & cath_proteins).length
			stats["COMBINED_RESULTS" + '_' + cath_tag + '_' + organism + "==funfamID_" + go_type +'--intersection_training_cath_FF (%)'] = (training_proteins & cath_proteins).length.fdiv(training_proteins.length)*100
			stats["COMBINED_RESULTS" + '_' + cath_tag + '_' + organism + "==funfamID_" + go_type +'--intersection_testing_cath_FF'] = (testing_proteins & cath_proteins).length
			stats["COMBINED_RESULTS" + '_' + cath_tag + '_' + organism + "==funfamID_" + go_type +'--intersection_testing_cath_FF (%)'] = (testing_proteins & cath_proteins).length.fdiv(testing_proteins.length)*100
		end
	end
	stats["COMBINED_RESULTS_" + organism + '==testing_proteins_with_geneID'] = testing_proteins.length
	stats["COMBINED_RESULTS_" + organism + '==testing_proteins_without_geneID (%)'] = testing_proteins_untranslated.length.fdiv(testing_proteins_untranslated.length + testing_proteins.length)*100
end

def statistics_report_data(container, html_file)
	template = File.open(File.join(REPORT_FOLDER, 'statistics_report.erb')).read
	report = Report_html.new(container, 'Full statistics for CAFA 2 and CAFA 3 datasets')
	report.build(template)
	report.write(html_file)
end

##########################
#OPT-PARSER
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:accessionID_dictionary] = nil
  opts.on("-d", "--accessionID_dictionary PATH", "Input file with gene accession ID - protein ID dictionary") do |data|
    options[:accessionID_dictionary] = data
  end

  options[:html_file] = "statistics_report.html"
  opts.on("-F", "--html_file PATH", "HTML file with statistics report for CAFA 2 and CAFA 3") do |html_file|
    options[:html_file] = html_file
  end

  
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
stats_complex = {}
geneAccession_protID_dictionary = {}
training_storage = {}
testing_storage = {}

load_geneAccession_protID_dictionary(options[:accessionID_dictionary], geneAccession_protID_dictionary)

tag_path_info.each do |id, path|
	if id.include?('input_cafa')
		get_input_cafa_data(path, id, geneAccession_protID_dictionary, stats_complex, training_storage, testing_storage)
	elsif id.include?('input_cath')
		get_input_CATH_data(path, id, stats_complex)
	elsif id.include?('input_networks')
		get_input_network_data(path, id, stats_complex)
	elsif id.include?('input_associations')
		info_string = 'assocs'
		get_input_assoc_predictions_data(path, id, info_string, stats_complex)
	elsif id.include?('input_predictions')
		info_string = 'preds'
		get_input_assoc_predictions_data(path, id, info_string, stats_complex)
	elsif id.include?('input_normpreds')
		info_string = 'norm_preds'
		get_input_assoc_predictions_data(path, id, info_string, stats_complex)
	end
end

#include this in get_combined_stats
# tag_path_info.each do |id, path|
# 	get_combined_stats(loaded_data, stats, path, geneAccession_protID_dictionary)			
# end

STDERR.puts JSON.pretty_generate(stats_complex)
Process.exit


File.open(options[:output_file], 'w') do |f|
	stats.each do |stat, values|
		f.puts "#{stat}\t#{values}"
	end
end

####- HTML REPORTING -####

# stats_ary = stats.to_a
# cafa2_human = stats_ary.select{|a| a.include?('CAFA2') && a.include?('HUMAN')}

# STDERR.puts cafa2_human.class
# Process.exit

# STDERR.puts JSON.pretty_generate(stats_complex)
# Process.exit

container = {
  :stats_info => stats.to_a
}

statistics_report_data(container, options[:html_file])

Process.exit