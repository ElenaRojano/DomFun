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

	testing_training_stats['training']['general_stats']['GO_terms'] = all_go_terms.uniq.length
	testing_training_stats['training']['general_stats']['uniq_proteins'] = annotated_proteins.map{|a| a.first}.uniq.length

	gomf_annotated_proteins = annotated_proteins.select{|a| a.last == 'F'}
	gobp_annotated_proteins = annotated_proteins.select{|a| a.last == 'P'}
	gocc_annotated_proteins = annotated_proteins.select{|a| a.last == 'C'}

	testing_training_stats['training']['GOMF']['proteins_annotated_in_GOMF'] = gomf_annotated_proteins.uniq.length
	testing_training_stats['training']['GOBP']['proteins_annotated_in_GOBP'] = gobp_annotated_proteins.uniq.length
	testing_training_stats['training']['GOCC']['proteins_annotated_in_GOCC'] = gocc_annotated_proteins.uniq.length
	

	testing_training_stats['training']['GOMF']['uniq_GOMF'] = go_subontologies['F'].uniq.length
	testing_training_stats['training']['GOBP']['uniq_GOBP'] = go_subontologies['P'].uniq.length
	testing_training_stats['training']['GOCC']['uniq_GOCC'] = go_subontologies['C'].uniq.length

	#TODO: protein_GOMF_relations and proteins_annotated_in give the same value.
	#Select which one to add.
	testing_training_stats['training']['GOMF']['protein_GOMF_relations'] = go_subontologies['F'].length
	testing_training_stats['training']['GOBP']['protein_GOBP_relations'] = go_subontologies['P'].length
	testing_training_stats['training']['GOCC']['protein_GOCC_relations'] = go_subontologies['C'].length

	return gomf_annotated_proteins, gobp_annotated_proteins, gocc_annotated_proteins
end

def calculate_testing_proteins_stats(testing_proteins, testing_training_stats, human_testing_training_stats, geneAccession_protID_dictionary, testing_storage)
	#TODO CONFIRMAR QUE LAS PROTEÍNAS DE TESTING EN EL HASH CORRESPONDEN A LA REALIDAD
	#HAY MAYOR NÚMERO DE PROTEÍNAS TRADUCIDAS QUE EN EL ARCHIVO ORIGINAL DE TESTING: ¿ES NORMAL?

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
	Dir.glob(path).each do |input_file|
		if input_file.include?('superfamilyID')
			tag = 'superfamilyID'
		elsif input_file.include?('funfamID')
			tag = 'funfamID'
		end		
		if input_file.include?('GOMF')
			annot = 'GOMF'
			#prefix = id.upcase + '==' + tag + '_' + annot + '_network--'
			#calculate_network_stats(input_file, stats, prefix, id, loaded_data, stats_complex, tag, annot)
			calculate_network_stats(input_file, id, stats_complex, tag, annot)
		elsif input_file.include?('GOBP')
			annot = 'GOBP'
			#prefix = id.upcase + '==' + tag + '_' + annot + '_network--'
			#calculate_network_stats(input_file, stats, prefix, id, loaded_data, stats_complex, tag, annot)
		elsif input_file.include?('GOCC')
			annot = 'GOCC'
			#prefix = id.upcase + '==' + tag + '_' + annot + '_network--'
			#calculate_network_stats(input_file, stats, prefix, id, loaded_data, stats_complex, tag, annot)
		end
	end

	#stats_complex['network'] = network_data

end

def calculate_network_stats(input_file, id, stats_complex, dom_tag, annot)
	network_data_complex = {'human' => {}, 'all' => {}}
	if input_file.include?('human')
		organism_tag = 'human'
	else
		organism_tag = 'all'
	end
	network_data_complex[organism_tag] = {dom_tag => {annot => {}}}
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
	#loaded_data['network_info'] = network_info
	#stats[prefix + 'uniq_proteins_number'] = network_info.keys.flatten.uniq.length
	#CHANGE TRAINING FOR NETWORK
	#CREATE NETWORK KEY BEFORE
	# STDERR.puts organism_tag.inspect
	# STDERR.puts stats_complex['network'][organism_tag].inspect
	# Process.exit
	network_data_complex[organism_tag][dom_tag][annot]['uniq_proteins_number'] = network_info.keys.flatten.uniq.length
	#stats[prefix + 'uniq_domains_number'] = network_info.values.flatten.select{|a| a.include?('.')}.uniq.length
	network_data_complex[organism_tag][dom_tag][annot]['uniq_domains_number'] = network_info.values.flatten.select{|a| a.include?('.')}.uniq.length
	#stats[prefix + 'uniq_GO_annotations_number'] = network_info.values.flatten.select{|a| a.include?('GO')}.uniq.length
	network_data_complex[organism_tag][dom_tag][annot]['uniq_GO_annotations_number'] = network_info.values.flatten.select{|a| a.include?('GO')}.uniq.length
	#stats[prefix + 'GO-proteins_relations'] = network_info.values.flatten.select {|el| el.include?('GO') }.length
	network_data_complex[organism_tag][dom_tag][annot]['GO-proteins_relations'] = network_info.values.flatten.select {|el| el.include?('GO') }.length
	#stats[prefix + 'domain-proteins_relations'] = network_info.values.flatten.select {|el| el.include?('.') }.length	
	network_data_complex[organism_tag][dom_tag][annot]['domain-proteins_relations'] = network_info.values.flatten.select {|el| el.include?('.') }.length
	#stats[prefix + 'total_domain/GO-protein_relations'] = network_info.values.flatten.length
	network_data_complex[organism_tag][dom_tag][annot]['total_domain/GO-protein_relations'] = network_info.values.flatten.length

	stats_complex['network'] = network_data_complex
	STDERR.puts JSON.pretty_generate(stats_complex)
	Process.exit
end



#Note: associations, predictions and normalized predictions
#load file are under the same method.
#Stats calculation has been joined in the same method aw.
def get_input_assoc_predictions_data(path, stats, id, info_string, loaded_data, stats_complex)
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
					calculate_stats(input_file, stats, prefix, id, mode, loaded_data)
				elsif input_file.include?('GOBP') && input_file.include?(assoc_method)
					prefix = id.upcase + '==superfamily_GOBP_' + assoc_method + info_string
 					calculate_stats(input_file, stats, prefix, id, mode, loaded_data)
				elsif input_file.include?('GOCC') && input_file.include?(assoc_method)
					prefix = id.upcase + '==superfamily_GOCC' + assoc_method + info_string
					calculate_stats(input_file, stats, prefix, id, mode, loaded_data)
				end
			end
		end	
	end
end

def calculate_stats(input_file, stats, prefix, id, mode, loaded_data)
	file_info = []
	open_file_save_array(input_file, file_info)
	if mode == 'associations'
		loaded_data['associations_info'] = file_info
		stats[prefix + 'GO-domain_relations_number'] = file_info.length
		stats[prefix + 'GO_number'] = file_info.map{|a| a.first}.uniq.length
		stats[prefix + 'domain_number'] = file_info.map{|a| a[1]}.uniq.length
		stats[prefix + 'max_assoc_val'] = file_info.map{|a| a.last}.max
		stats[prefix + 'min_assoc_val'] = file_info.map{|a| a.last}.min
		assoc_vals = file_info.map{|a| a.last.to_f}
		stats[prefix + 'ave_assoc_val'] = assoc_vals.inject { |sum, n| sum + n }.fdiv(assoc_vals.length)
	elsif mode == 'predictions'
		loaded_data['predictions_info'] = file_info
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
		loaded_data['normalized_predictions_info'] = file_info
		stats[prefix + 'predicted_proteins_after_norm'] = file_info.map{|a| a.first}.uniq.length
		stats[prefix + 'total_GOs_predicted_after_norm'] = file_info.map{|a| a[2]}.uniq.length
		stats[prefix + 'total_domains_after_norm'] = file_info.map{|a| a[1].split(',')}.flatten.uniq.length
		stats[prefix + 'max_norm_pred_val'] = file_info.map{|a| a.last}.max
		stats[prefix + 'min_norm_pred_val'] = file_info.map{|a| a.last}.min
		pred_vals = file_info.map{|a| a.last.to_f}
		stats[prefix + 'ave_norm_pred_val'] = pred_vals.inject { |sum, n| sum + n }.fdiv(pred_vals.length)
	end
end

def get_combined_stats(loaded_data, stats, path, geneAccession_protID_dictionary)
	#loaded_data['associations_info'] = asociaciones
	#loaded_data['predictions_info']	= predicciones
	#loaded_data['normalized_predictions_info'] = predicciones normalizadas
	#loaded_data['network_info'] = informacion_red
	#loaded_data['training_proteins'] = CAFA training
	#loaded_data['training_proteins_human'] = CAFA training_human
	#loaded_data['testing_proteins'] = CAFA testing (in gene accession ID)
	#loaded_data['testing_proteins_translated'] = testing proteins in proteinID
	#loaded_data['cath_info_SF'] = CATH info (SF)
	#loaded_data['cath_info_FF'] = CATH info (FF)
	#loaded_data['testing_proteins_translated'] = testing proteins translated to UniProt ID
	#loaded_data['testing_proteins_untranslated'] = testing proteins without translation to UniProt ID

	#geneAccession_protID_dictionary == {P32234 =>128UP_DROME}
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
stats = {}
stats_complex = {}
loaded_data ={}
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
		info_string = '_assocs--'
		get_input_assoc_predictions_data(path, stats, id, info_string, loaded_data, stats_complex)
	elsif id.include?('input_predictions')
		info_string = '_preds--'
		get_input_assoc_predictions_data(path, stats, id, info_string, loaded_data, stats_complex)
	elsif id.include?('input_normpreds')
		info_string = info_string = '_norm_preds--'
		get_input_assoc_predictions_data(path, stats, id, info_string, loaded_data, stats_complex)
	end
end

#include this in get_combined_stats
tag_path_info.each do |id, path|
	get_combined_stats(loaded_data, stats, path, geneAccession_protID_dictionary)			
end


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

STDERR.puts JSON.pretty_generate(stats_complex)
Process.exit

container = {
  :stats_info => stats.to_a
}

statistics_report_data(container, options[:html_file])

Process.exit