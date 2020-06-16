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
require 'optparse'
require "statistics2"
require 'bigdecimal'


##########################
#METHODS
##########################

def load_predictions_file(predictions_file)
	predictions = []
	File.open(predictions_file).each do |line|
		line.chomp!
		next if line.include?('ProteinID')
		protein, domains, funSys, combinedScore = line.split("\t")
		predictions << [protein, funSys, combinedScore.to_f]
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

def load_prediction(pairs_array)
	pred = {}
	min = nil
	max = nil
	pairs_array.each do |key, label, score| #protein, FunSys, assocScore
		query = pred[key]
		if !min.nil? && !max.nil?
			min = score if score < min
			max = score if score > max
		else
			min = score; max = score
		end
		if query.nil?
			pred[key] = [[label], [score]]
		else
			query.first << label
			query.last << score
		end
	end
	return pred, [min, max]
end


# Pandey 2007, Association Analysis-based Transformations for Protein Interaction Networks: A Function Prediction Case Study
def get_pred_rec(meth, cut_number = 100, top_number = 10000, control_protein_FunSys, predictions)
	performance = [] #cut, pred, rec
	preds, limits = load_prediction(predictions)
	cuts = get_cuts(limits, cut_number)
	cuts.each do |cut|
		prec, rec = pred_rec(preds, cut, top_number, control_protein_FunSys)
		performance << [cut, prec, rec]
	end
	return performance
end

def pred_rec(preds, cut, top, control_protein_FunSys)
	predicted_labels = 0 #m
	true_labels = 0 #n
	common_labels = 0 # k
	control_protein_FunSys.each do |key, c_labels|
		true_labels += c_labels.length #n
		pred_info = preds[key]
		if !pred_info.nil?
			labels, scores = pred_info
			reliable_labels = get_reliable_labels(labels, scores, cut, top)
			predicted_labels += reliable_labels.length #m
			common_labels += (c_labels & reliable_labels).length #k
		end
	end
	#puts "cut: #{cut} trueL: #{true_labels} predL: #{predicted_labels} commL: #{common_labels}"
	prec = common_labels.to_f/predicted_labels
	rec = common_labels.to_f/true_labels
	prec = 0.0 if prec.nan?
	rec = 0.0 if rec.nan?
	return prec, rec
end

def get_cuts(limits, n_cuts)
	cuts = []
	range = (limits.last - limits.first).abs.fdiv(n_cuts)
	range = BigDecimal(range, 10)
	cut = limits.first
	(n_cuts + 1).times do |n|
		cuts << (cut + n * range).to_f
	end
	return cuts
end

def get_reliable_labels(labels, scores, cut, top)
	reliable_labels = []
	scores.each_with_index do |score, i|
		reliable_labels << [labels[i], score] if score >= cut
	end
	reliable_labels = reliable_labels.sort!{|l1,l2| l2.last <=> l1.last}[0..top-1].map{|pred| pred.first}
	return reliable_labels
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

	options[:assoc_meth] = nil
		opts.on("-m", "--assoc_meth STRING", "Association method used") do |data|
		options[:assoc_meth] = data
	end
	
	options[:performance_file] = 'precision_recall.txt'
	opts.on("-p", "--performance_file PATH", "Output file with PR values") do |data|
	options[:performance_file] = data
	end

	opts.on_tail("-h", "--help", "Show this message") do
		puts opts
		exit
	end

end.parse!

##########################
#MAIN
##########################

control_protein_FunSys = load_control_file(options[:control_file])

domains_FunSys_predictions = load_predictions_file(options[:input_predictions])

performance = get_pred_rec(options[:assoc_meth], 100, 10000, control_protein_FunSys, domains_FunSys_predictions)

File.open(options[:performance_file], 'w') do |f|
	f.puts %w[cut prec rec meth].join("\t")
	performance.each do |item|
		item << options[:assoc_meth].to_s
		f.puts item.join("\t")
	end
end