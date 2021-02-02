#! /usr/bin/env ruby

##########################
# Rojano E. & Seoane P., Feb 2021
# Generate files for CAFA 3 validation (from normalized predictions)
##########################

REPORT_FOLDER=File.expand_path(File.join(File.dirname(__FILE__), '..', 'templates'))
ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'DomFun'))
require 'generalMethods'
require 'optparse'
#require 'csv'

##########################
#METHODS
#########################
# Corregir el input de combined scores

def load_predictions(input_file)
	predictions = {}
	File.open(input_file).each do |line|
		line.chomp! 
    proteinID, domains, go_term, p_value = line.split("\t")	
    query = predictions[proteinID]
    if query.nil?
      predictions[proteinID] = [[go_term, p_value.to_f]]
    else
      query << [go_term, p_value.to_f]
    end 
  end
  return predictions
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

##########################
#OPT-PARSER
##########################
options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:input_predictions] = nil
  opts.on("-a", "--input_predictions PATH", "Input file with normalized predictions") do |input_predictions|
    options[:input_predictions] = input_predictions
  end

  options[:cath_dict] = nil
  opts.on("-c", "--cath_dict PATH", "Input CATH file to translate UniProtIDs to GeneIDs") do |data|
    options[:cath_dict] = data
  end

  options[:cafa_dict] = nil
  opts.on("-d", "--cafa_dict PATH", "Input CAFA3 file to translate GeneIDs to CAFA3IDs") do |data|
    options[:cafa_dict] = data
  end  

  options[:domains_class] = 'funfamID'
  opts.on("-D", "--domains_class STRING", "Domains type") do |data|
    options[:domains_class] = data
  end  

  options[:output_file] = 'results_to_CAFA3_validation.txt'
  opts.on("-o", "--output_file PATH", "Output file") do |output_file|
    options[:output_file] = output_file
  end

  options[:untranslated_proteins] = 'untranslated_proteins.txt'
  opts.on("-u", "--untranslated_proteins PATH", "Output file with untranslated proteins list") do |data|
    options[:untranslated_proteins] = data
  end

end.parse!

##########################
#MAIN
##########################

predictions = load_predictions(options[:input_predictions])
cath_dict = load_hash(options[:cath_dict], 'a')
cafa_dict = load_hash(options[:cafa_dict], 'b')
cafa3_predictions, untranslated_proteins = translate_uniprot_to_CAFA(predictions, cath_dict, cafa_dict)

handler = File.open(options[:output_file], 'w')
cafa3_predictions.each do |cafaID, predictions|
	predictions.each do |prediction|
		handler.puts "#{cafaID}\t#{prediction.join("\t")}"
	end
end

handler = File.open(options[:untranslated_proteins], 'w')
handler.puts untranslated_proteins.join("\n")