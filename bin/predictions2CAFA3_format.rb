#! /usr/bin/env ruby

##########################
# Rojano E. & Seoane P., Feb 2021
# Generate files for CAFA 3 validation (from normalized predictions)
# Load several files with normalized predictions and for each one returns:
# 1. A yaml file for each predictions file included with info for CAFA 3 validation
# 2. A yaml file with info for all predictions to plot CAFA 3 validation graphs
##########################

REPORT_FOLDER=File.expand_path(File.join(File.dirname(__FILE__), '..', 'templates'))
ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'DomFun'))
require 'generalMethods'
require 'optparse'
require 'yaml'
require 'fileutils'
#require 'csv'

##########################
#METHODS
#########################
# Corregir el input de combined scores

def load_predictions(prediction_files)
	all_predictions = {}
  predictions = {}
	prediction_files.each do |f|
    filename = File.basename(f).split(".").first
    File.open(f).each do |line|
  		line.chomp! 
      proteinID, domains, go_term, p_value = line.split("\t")	
      p_value = sprintf "%.2f", p_value.to_f
      query = predictions[proteinID]
      if query.nil?
        predictions[proteinID] = [[go_term, p_value]]
      else
        query << [go_term, p_value]
      end 
    end
    all_predictions[filename] = predictions #SAVE ALL HASHES FOR EACH PREDICTION FILE
  end
  return all_predictions
end

def translate_uniprot_to_CAFA(predictions, cath_dict, cafa_dict)
	cafa3_predictions = {}
	untranslated_proteins = []
	predictions.each do |pred_UniprotID, content|
		if cath_dict[pred_UniprotID].nil?
			untranslated_proteins << pred_UniprotID
		else
			geneID = cath_dict[pred_UniprotID].first
		end
		if cafa_dict[geneID].nil?
			untranslated_proteins << pred_UniprotID
		else
			cafaID = cafa_dict[geneID].first
			cafa3_predictions[cafaID] = content
		end
	end
	return cafa3_predictions, untranslated_proteins
end

def generate_yaml_hash_predictions(output_file_path, obo_path, benchmark_files, yaml_output)
  yaml_hash = {
    "assess" => {
      "file" => output_file_path, "obo" => obo_path, "benchmark" => benchmark_files, "results" => yaml_output
    }
  }
  return yaml_hash
end

def generate_plot_yaml(plot_results_path, title, smooth, prediction_files)
  filenames = prediction_files.map do |f|
    File.basename(f).split(".").first
  end
  yaml_hash = {
    "plot" => {
      "results" => plot_results_path, "title" => title, "smooth" => smooth
    }
  }
  counter = 1
  filenames.each do |filename|
    yaml_hash["plot"]["file" + counter.to_s] = filename
    counter += 1
  end
  return yaml_hash

end

def generate_output_path(domain_class, association_method, string, file_path)
  output_name = []
  output_name << domain_class
  output_name << association_method
  output_name << string
  filename = output_name.join('_')
  full_output_path = File.join(file_path, filename)
  return full_output_path, filename
end

# def parse_output_file()

# end

##########################
#OPT-PARSER
##########################
options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:input_predictions] = nil
  opts.on("-a", "--input_predictions PATH", "Input files with normalized predictions") do |input_predictions|
    options[:input_predictions] = input_predictions
  end

  options[:benchmark_path] = nil
  opts.on("-b", "--benchmark_path PATH", "Input PATH with benchmark files for CAFA 3") do |data|
    options[:benchmark_path] = data
  end  

  options[:accesion_geneid_dictionary] = nil
  opts.on("-c", "--accesion_geneid_dictionary PATH", "Input UniProt accesion-GeneID dictionary file to translate UniProtIDs to GeneIDs") do |data|
    options[:accesion_geneid_dictionary] = data
  end

  options[:targetID_geneid_dictionary] = nil
  opts.on("-d", "--targetID_geneid_dictionary PATH", "Input file to translate GeneIDs to CAFA3 IDs") do |data|
    options[:targetID_geneid_dictionary] = data
  end  

  options[:domains_class] = 'funfamID'
  opts.on("-e", "--domains_class STRING", "Domains type to use as group name in yaml file") do |data|
    options[:domains_class] = data
  end  

  options[:association_method] = 'hypergeometric'
  opts.on("-f", "--association_method STRING", "Association method to use as mode in yaml file") do |data|
    options[:association_method] = data
  end  

  options[:cafa3_results] = 'results4CAFA'
  opts.on("-g", "--cafa3_results PATH", "Path to save CAFA 3 plot files") do |data|
    options[:cafa3_results] = data
  end

  options[:set_model] = '1'
  opts.on("-m", "--set_model STRING", "Set model for CAFA 3 analysis") do |data|
    options[:set_model] = data
  end

  options[:output_path] = nil
  opts.on("-o", "--output_path PATH", "Output file") do |data|
    options[:output_path] = data
  end

  options[:obo_path] = nil
  opts.on("-p", "--obo_path PATH", "OBO file for CAFA 3 analysis") do |data|
    options[:obo_path] = data
  end

  options[:untranslated_proteins] = nil
  opts.on("-u", "--untranslated_proteins PATH", "Output path with untranslated proteins list") do |data|
    options[:untranslated_proteins] = data
  end

  options[:yaml_file] = nil
  opts.on("-y", "--yaml_file PATH", "Output path to storage yaml files for CAFA 3 validation") do |data|
    options[:yaml_file] = data
  end

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!

##########################
#MAIN
##########################
prediction_files = Dir.glob(options[:input_predictions])
all_predictions = load_predictions(prediction_files)
accesion_geneid_dictionary = load_hash(options[:accesion_geneid_dictionary], 'a')
targetID_geneid_dictionary = load_hash(options[:targetID_geneid_dictionary], 'b')

# OUTPUT PREDICTION FILE FORMAT FOOTER:
footer = "END"

# YAML FILES PARAMETERS:
obo_path = options[:obo_path]
benchmark_files = options[:benchmark_path]
title = "all"
smooth = "N"

domain_class = nil
association_method = nil
FileUtils.mkdir_p(options[:cafa3_results])

all_prediction_filenames = []
# REFORMED CAFA 3 PREDICTION FILES:
all_predictions.each do |filename, prediction|
  cafa3_predictions, untranslated_proteins = translate_uniprot_to_CAFA(prediction, accesion_geneid_dictionary, targetID_geneid_dictionary)
  domain_class = filename.split('_')[1]
  association_method = filename.split('_')[2]
  header = "AUTHOR #{domain_class + association_method}\nMODEL #{options[:set_model]}\nKEYWORDS sequence alignment.\n"
  cafa_name = domain_class + association_method
  full_output_path, filename = generate_output_path(cafa_name, options[:set_model], 'all.txt', options[:output_path])
  all_prediction_filenames << filename

  File.open(full_output_path, 'w') do |f|
    f.puts header
    cafa3_predictions.each do |cafaID, predictions|
    	predictions.each do |prediction|
    		f.puts "#{cafaID}\t#{prediction.join("\t")}"
    	end
    end
    f.puts footer
  end

  # GENERATE YAML FILES FOR PREDICTIONS:
  
  output_file_path = full_output_path #where CAFA 3 predictions file is stored
  full_output_path, filename = generate_output_path(domain_class, association_method, 'config_launch.yaml', options[:yaml_file])
  output_file1 = [options[:domains_class], options[:set_model], title].join('_')


  yaml_pred_hash = generate_yaml_hash_predictions(output_file_path, obo_path, benchmark_files, options[:cafa3_results])
  File.open(full_output_path, "w") { |file| file.write(yaml_pred_hash.to_yaml) }

  # GENERATE FILE WITH LOST PROTEINS:
  full_output_path, filename = generate_output_path(domain_class, association_method, 'unnanotated_proteins.txt', options[:untranslated_proteins])
  File.open(full_output_path, 'w') do |f|
    f.puts untranslated_proteins.join("\n")
  end

end

#Generate uniq PLOT file
yaml_plot_hash = generate_plot_yaml(options[:cafa3_results], title, smooth, all_prediction_filenames)
yaml_plot_all_filename = File.join(options[:yaml_file], 'plot_all.yaml')
File.open(yaml_plot_all_filename, "w") { |file| file.write(yaml_plot_hash.to_yaml) }