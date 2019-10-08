#! /usr/bin/env ruby

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

def normalize_association_values(cafa_predictions)
  raw_association_values = []
  normalized_prediction_values = []
  cafa_predictions.each do |protID, goTerm, value|
    raw_association_values << value.to_f
  end
  raw_values_mean = raw_association_values.mean
  raw_values_sd = raw_association_values.standard_deviation
  association_values_minus_mean = []
  raw_association_values.each do |value|
    z_score = (value - raw_values_mean).fdiv(raw_values_sd)
    if z_score > 2 
      z_score = 2
    elsif z_score < -2
      z_score = -2
    end
    z_score = z_score.fdiv(5) + 0.5
    association_values_minus_mean << z_score
  end
  cafa_predictions.each_with_index do |info, c|
    protID, goTerm, association_values = info
    normalized_prediction_values << [protID, goTerm, association_values_minus_mean[c]]
  end
  cafa_predictions = normalized_prediction_values
  return cafa_predictions
end

##########################
#OPT-PARSER
##########################
options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:input_predictions] = nil
  opts.on("-a", "--input_predictions PATH", "Input predictions file") do |input_predictions|
    options[:input_predictions] = input_predictions
  end

  options[:input_cafa] = nil
  opts.on("-c", "--input_cafa PATH", "Input CAFA file to translate UniProtIDs to CAFAIDs") do |input_cafa|
    options[:input_cafa] = input_cafa
  end

  options[:output_file] = 'results_to_CAFA2_validation.txt'
  opts.on("-o", "--output_file PATH", "Output file") do |output_file|
    options[:output_file] = output_file
  end

  options[:do_norm] = false
  opts.on("-n", "--do_norm", "Normalize prediction values") do
    options[:do_norm] = true
  end

end.parse!

##########################
#MAIN
##########################

predictions = load_predictions(options[:input_predictions])
cafa_predictions = translate_protein_ids(options[:input_cafa], predictions)
cafa_predictions = normalize_association_values(cafa_predictions) if options[:do_norm]

handler = File.open(options[:output_file], 'w')
cafa_predictions.each do |info|
  handler.puts "#{info.join("\t")}"
end