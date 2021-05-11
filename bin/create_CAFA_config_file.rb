#! /usr/bin/env ruby

##########################
# Rojano E. & Seoane P., July 2019
# Script to generate the config file for CAFA assessment
##########################

require 'optparse'

##########################
#METHODS
##########################

##########################
#OPT-PARSER
##########################
options = {}
OptionParser.new do |opts|
	opts.banner = "Usage: #{__FILE__} [options]"

	options[:assessment_dir] = nil
		opts.on("-a", "--assessment_dir PATH", "Assessment results directory") do |data|
		options[:assessment_dir] = data
	end

	options[:benchmark_type] = '2'
		opts.on("-b", "--benchmark_type STRING", "Evaluation benchmark type. 1: no-knowledge (NK); 2: limited-knowledge (LK)") do |data|
		options[:benchmark_type] = data
	end

	options[:benchmark] = 'benchmark'
		opts.on("-B", "--benchmark STRING", "Benchmark folder name") do |data|
		options[:benchmark] = data
	end

	options[:cafa_path] = nil
		opts.on("-c", "--cafa_path PATH", "CAFA files (benchmark groundtruth and lists path") do |data|
		options[:cafa_path] = data
	end

	options[:evaluation_mode] = '2'
		opts.on("-e", "--evaluation_mode STRING", "Evaluation mode. 1: full; 2: partial") do |data|
		options[:evaluation_mode] = data
	end

	options[:ontology_folder] = 'ontology'
		opts.on("-f", "--ontology_folder STRING", "Ontology folder name") do |data|
		options[:ontology_folder] = data
	end

	options[:go_subontology] = ['mfo']
		opts.on("-g", "--go_subontology STRING", "GO subontologies: mfo, bpo or cco") do |data|
		options[:go_subontology] = data.split(',')
	end

	options[:config_output] = 'config.job'
		opts.on("-o", "--config_output PATH", "Config job output file") do |data|
		options[:config_output] = data
	end

	options[:performance_metrics] = ['f']
		opts.on("-p", "--performance_metrics STRING", "Evaluation performance. Accepts multiple selection (comma separated). f (F1-max), s (S2-min), wf (weighted F1-max), ns (normalized S2-min) and auc (area under the ROC curve).") do |data|
		options[:performance_metrics] = data.split(',')
	end

	options[:register_filename] = 'register.tab'
		opts.on("-r", "--register_filename STRING", "Register filename") do |data|
		options[:register_filename] = data
	end

	options[:sample_asignation] = 'files_to_cafa'
		opts.on("-s", "--sample_asignation STRING", "Sample assignation filename") do |data|
		options[:sample_asignation] = data
	end

	options[:taxon_category] = 'all'
		opts.on("-t", "--taxon_category STRING", "Organism to assess predictions. Please choose between all, easy, hard, eukarya, prokarya, or SPECIES (e.g., HUMAN, MOUSE)") do |data|
		options[:taxon_category] = data
	end

	options[:general_variables] = 'general_variables.txt'
		opts.on("-V", "--general_variables PATH", "File with general variables for CAFA analysis") do |data|
		options[:general_variables] = data
	end

	opts.on_tail("-h", "--help", "Tool information") do
		puts opts
		exit
	end

end.parse!

##########################
#MAIN
##########################

if options[:benchmark_type] == '1'
	label_benchmark_type = 'type1'
	#label_benchmark_type = 'NK'
elsif options[:benchmark_type] == '2'
	label_benchmark_type = 'type2'
	#label_benchmark_type = 'LK'
else
	abort('Wrong argument. Please select between 1: no-knowledge (NK); 2: limited-knowledge (LK)')
end

config_filenames = []
options[:go_subontology].each do |subontology|
	lists_filename = [subontology, options[:taxon_category], label_benchmark_type].join('_') + '.txt'
	bootstrap = [subontology, options[:taxon_category], 'type' + options[:benchmark_type]].join('_') + '.mat'
	filepath = options[:config_output] + subontology + '.job'
	config_filenames << File.basename(filepath)
	File.open(filepath, 'w') do |f|
		f.puts "pred_dir = #{File.join(options[:assessment_dir], 'prediction', subontology)}" 
		f.puts "prev_dir = #{File.join(options[:assessment_dir], 'seq-centric', subontology)}" 
		f.puts "eval_dir = #{File.join(options[:assessment_dir], 'evaluation')}" 
		f.puts "annotation = #{File.join(options[:cafa_path], options[:benchmark], 'groundtruth', subontology + 'a.mat')}" 
		f.puts "benchmark = #{File.join(options[:cafa_path], options[:benchmark], 'lists', lists_filename)}" 
		f.puts "bootstrap = #{File.join(options[:assessment_dir], 'bootstrap', bootstrap)}" 
		f.puts "ontology = #{subontology}"
		f.puts "category = #{options[:taxon_category]}"
		f.puts "type = #{options[:benchmark_type]}"
		f.puts "mode = #{options[:evaluation_mode]}"
		options[:performance_metrics].each do |metric|
			f.puts "metric = #{metric}"
		end
		f.puts "model = all"
		f.puts "model = +BN4S"
		f.puts "model = +BB4S"
	end
end


ontologies = options[:go_subontology].map{|o| o.upcase}.join(',')
config_files = config_filenames.join(',')
File.open(File.join(options[:cafa_path], 'matlab', options[:general_variables]), 'w') do |f|
	f.puts "cafa_system_dir = '#{options[:cafa_path]}'" 
	f.puts "dir_validation = '#{options[:assessment_dir]}'"
	f.puts "ontology_folder = '#{File.join(options[:cafa_path], options[:ontology_folder])}'"
	f.puts "ontologies = '#{ontologies}'"
	f.puts "dir_evaluation_configs = '#{File.dirname(options[:config_output])}'"
	f.puts "config_files = '#{config_files}'"
	f.puts "register_file = '#{options[:register_filename]}'"
	f.puts "sample_asignation = '#{options[:sample_asignation]}'"
end
