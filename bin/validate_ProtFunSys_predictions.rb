#! /usr/bin/env ruby
##########################
# Rojano E. & Seoane P., June 2019
# Domain to functional annotation predictor validation system
# The script uses the predictions file and proteins-FunSys from UniProtKB
# It compares the predictions with the proteins-FunSys to validate the functioning of the predictor
# Generate values to plot in a PR
##########################


REPORT_FOLDER=File.expand_path(File.join(File.dirname(__FILE__), '..', 'templates'))
ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'DomFun'))
require 'generalMethods.rb'
require 'csv'
require 'optparse'
require "statistics2"
require "terminal-table"
require 'report_html'


##########################
#METHODS
##########################

def load_predictions_file(predictions_file)
	predictions = {}
	File.open(predictions_file).each do |line|
		line.chomp!
		next if line.include?('ProteinID')
		protein, domains, funSys, assocVal, combinedScore = line.split("\t")
		query = predictions[protein]
		if query.nil?
			predictions[protein] = [funSys]
		else
			query << funSys
		end
	end
	return predictions
end

def load_control_file(control_file)
	control_protein_FunSys = {}
	File.open(control_file).each do |line|
		line.chomp!
		proteinID, funSys = line.split("\t")
		control_protein_FunSys[proteinID] = funSys.split(';')
	end
	return control_protein_FunSys
end

##########################
#OPT-PARSER
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"
  
	options[:input_predictions] = nil
	opts.on("-a", "--input_predictions PATH", "Domain-function predictions") do |data|
	options[:input_predictions] = data
	end

	options[:control_file] = nil
	opts.on("-c", "--control_file PATH", "Control dataset with proteins-FunSys from UniProtKB") do |data|
	options[:control_file] = data
	end

	opts.on_tail("-h", "--help", "Show this message") do
	puts opts
	exit
  end

end.parse!

##########################
#MAIN
##########################

domains_FunSys_predictions = load_predictions_file(options[:input_predictions])
control_protein_FunSys = load_control_file(options[:control_file])

# domains_FunSys_predictions.keys.inspect
number_of_predicted_FunSys = 0
number_of_FunSys_per_protein = 0
number_of_shared_FunSys = 0
control_protein_FunSys.each do |protein, funSys|
	predicted_FunSys = domains_FunSys_predictions[protein]
	number_of_FunSys_per_protein += funSys.length
	unless predicted_FunSys.nil?
		number_of_predicted_FunSys += predicted_FunSys.length
		if !predicted_FunSys.nil?
			commonFunsys = predicted_FunSys & funSys
			number_of_shared_FunSys += commonFunsys.length
		end
	end
end
precision = number_of_shared_FunSys.fdiv(number_of_predicted_FunSys)
STDERR.puts precision
recall = number_of_shared_FunSys.fdiv(number_of_FunSys_per_protein)
STDERR.puts recall
