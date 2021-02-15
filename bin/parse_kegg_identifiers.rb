#! /usr/bin/env ruby

##########################
# Rojano E. & Seoane P., Feb 2021
# Parse KEGG IDs
# KEGG ontology points at the same pathway using different identifiers (depending on the specie)
# With this script we ensure we have the same pathway id for all species
# For Reactome, we directly translate IDs in add_protein_functional_families.rb, before printing tripartite networks
##########################

REPORT_FOLDER=File.expand_path(File.join(File.dirname(__FILE__), '..', 'templates'))
ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'DomFun'))
require 'generalMethods'
require 'optparse'

##########################
#METHODS
#########################

def translate_terms(kegg_data)
	kegg_data.transform_values!{|v| v.map{|id| id.gsub(/path:([a-z]{3,4})/, 'path:map') unless id.nil?}}
end

##########################
#OPT-PARSER
##########################
options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:kegg_file] = nil
  opts.on("-a", "--kegg_file PATH", "Input KEGG file to translate IDs") do |data|
    options[:kegg_file] = data
  end

  options[:output_kegg] = 'output_kegg.txt'
  opts.on("-o", "--output_kegg PATH", "Output file with KEGG IDs translated") do |data|
    options[:output_kegg] = data
  end

end.parse!

##########################
#MAIN
##########################

kegg_data = load_hash(options[:kegg_file], 'a')
translate_terms(kegg_data)
File.open(options[:output_kegg], 'w') do |f|
	kegg_data.each do |gene, pathways|
		pathways.each do |pathway|
			f.puts "#{gene}\t#{pathway}"
		end
	end
end