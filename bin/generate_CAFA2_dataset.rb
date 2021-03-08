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

##########################
#MODULES
#########################

module Enumerable

    def sum
      self.inject(0){|accum, i| accum + i }
    end

    def mean
      self.sum/self.length.to_f
    end

    def sample_variance
      m = self.mean
      sum = self.inject(0){|accum, i| accum +(i-m)**2 }
      sum/(self.length - 1).to_f
    end

    def standard_deviation
      Math.sqrt(self.sample_variance)
    end

end 

##########################
#METHODS
#########################
# Corregir el input de combined scores

def load_predictions(input_files)
	predictions = {}
  Dir.glob(input_files).each do |input_file|
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
  end
  return predictions
end

def translate_protein_ids(cafa_file, predictions)
  cafa_predictions = []
  File.open(cafa_file).each do |line|
    line.chomp!
    cafa_data = line.split("\t")
    cafaID = cafa_data[3]
    proteinID = cafa_data[6]
    go_associations = predictions[proteinID]
    unless go_associations.nil?
      go_associations.each do |goID, value|
        cafa_predictions << [cafaID, goID, value]
      end
    end 
  end
  return cafa_predictions
end

def translate_uniprot_to_CAFA3(predictions, accesion_geneid_dictionary, targetID_geneid_dictionary)
  #{"U3KPV4"=>[["GO:0030259", 0.999999755632719], ["GO:0009624", 0.9922574199563998]]}
  cafa3_predictions = {}
  predictions.each do |pred_UniprotID, content|
    if !accesion_geneid_dictionary[pred_UniprotID].nil?
      geneID = accesion_geneid_dictionary[pred_UniprotID].first
      if !targetID_geneid_dictionary[geneID].nil?
        cafaID = targetID_geneid_dictionary[geneID].first
        cafa3_predictions[cafaID] = content
      end
    end
  end
  return cafa3_predictions
end


##########################
#OPT-PARSER
##########################
options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:input_predictions] = nil
  opts.on("-a", "--input_predictions PATH", "Input predictions file") do |data|
    options[:input_predictions] = data
  end

  options[:input_dictionary] = nil
  opts.on("-c", "--input_dictionary PATH", "Input CAFA file to translate UniProtIDs to CAFAIDs (For CAFA 2)") do |data|
    options[:input_dictionary] = data
  end

  options[:accesion_geneid_dictionary] = nil
  opts.on("-g", "--accesion_geneid_dictionary PATH", "Uniprot accesion - Geneid dictionary") do |data|
    options[:accesion_geneid_dictionary] = data
  end

  options[:cafa_mode] = 'CAFA2'
  opts.on("-m", "--cafa_mode STRING", "CAFA mode to generate file for assessment") do |data|
    options[:cafa_mode] = data
  end

  options[:output_file] = 'results_to_CAFA2_validation.txt'
  opts.on("-o", "--output_file PATH", "Output file") do |output_file|
    options[:output_file] = output_file
  end

  options[:targetID_geneid_dictionary] = nil
  opts.on("-t", "--targetID_geneid_dictionary PATH", "Target ID - Geneid dictionary") do |data|
    options[:targetID_geneid_dictionary] = data
  end

end.parse!

##########################
#MAIN
##########################

predictions = load_predictions(options[:input_predictions])
if options[:cafa_mode] == 'CAFA2'
  cafa_predictions = translate_protein_ids(options[:input_dictionary], predictions)
  File.open(options[:output_file], 'w') do |f|
    cafa_predictions.each do |info|
      f.puts "#{info.join("\t")}"
    end
  end
elsif options[:cafa_mode] == 'CAFA3'
  accesion_geneid_dictionary = load_hash(options[:accesion_geneid_dictionary], 'a')
  targetID_geneid_dictionary = load_hash(options[:targetID_geneid_dictionary], 'b')
  cafa_predictions = translate_uniprot_to_CAFA3(predictions, accesion_geneid_dictionary, targetID_geneid_dictionary)
  File.open(options[:output_file], 'w') do |f|
    cafa_predictions.each do |targetID, info|
      info.each do |i|
        f.puts "#{targetID}\t#{i.join("\t")}"
      end
    end
  end
else
  abort('Wrong CAFA mode. Please choose between CAFA2 or CAFA3')
end