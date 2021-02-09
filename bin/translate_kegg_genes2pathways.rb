#! /usr/bin/env ruby
##########################
# Rojano E. & Seoane P., July 2019
# Translate KEGG genes into pathways
##########################

REPORT_FOLDER=File.expand_path(File.join(File.dirname(__FILE__), '..', 'templates'))
ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'DomFun'))
# require 'generalMethods.rb'
require 'optparse'

##########################
#METHODS
##########################

def load_kegg_dictionary(pathway_to_genes_file)
	kegg_dictionary = {}
	File.open(pathway_to_genes_file).each do |line|
		line.chomp!
		keggGeneID, keggPathwayID = line.split("\t")
		query = kegg_dictionary[keggGeneID]
		if query.nil?
			kegg_dictionary[keggGeneID] = [keggPathwayID]
		else
			query << keggPathwayID
		end

	end
	return kegg_dictionary
end

def load_network(network_kegg_file)
	network_kegg = []
	domain_ids = []
	File.open(network_kegg_file).each do |line|
		line.chomp!
		item, gene = line.split("\t")
		if !item.match(/^[a-z]{3,3}:/).nil?
			network_kegg << [item, gene]
		else
			domain_ids << [item, gene]
		end
	end
	return network_kegg, domain_ids
end

##########################
#OPT-PARSER
##########################
options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:kegg_pathways] = nil
  opts.on("-k", "--kegg_pathways PATH", "Input file with KEGG genes to pathways") do |data|
    options[:kegg_pathways] = data
  end

  options[:network_kegg] = nil
  opts.on("-n", "--network_kegg PATH", "Network with KEGG genes to translate into pathways") do |data|
    options[:network_kegg] = data
  end

  options[:output_path] = 'network_kegg_pathways'
  opts.on("-o", "--output_path PATH", "Resulting network output path") do |data|
    options[:output_path] = data
  end

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!

##########################
#MAIN
##########################
kegg_dictionary = load_kegg_dictionary(options[:kegg_pathways])
network_kegg, domain_ids = load_network(options[:network_kegg])
pathways_network = []
network_kegg.each do |kegg, proteinID|
	pathwayIDs = kegg_dictionary[kegg]
	unless pathwayIDs.nil?
		pathwayIDs.each do |pathway|
			pathways_network << [pathway, proteinID]
		end
	end
end
File.open(options[:output_path], 'w') do |f|
	pathways_network.each do |pair|
		f.puts pair.join("\t")
	end
	domain_ids.each do |pair|
		f.puts pair.join("\t")
	end
end