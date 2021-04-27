#! /usr/bin/env ruby
##########################
# Rojano E. & Seoane P., July 2019
# Generate training and testing datasets from CAFA2 and UniProt data.
# This script translate identifiers and generate the file for the
# construction of tripartite network (domain-protein-FunSys)
##########################

REPORT_FOLDER=File.expand_path(File.join(File.dirname(__FILE__), '..', 'templates'))
ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'DomFun'))
require 'generalMethods.rb'
require 'optparse'
require 'csv'
require 'semtools'

##########################
#METHODS
##########################

def build_tripartite_network(domain_tuples, annot_tuples, filename)
	handler = File.open(filename, 'w')
	annot_tuples.each do |protID, annots|
		annots.each do |a|
			handler.puts "#{a}\t#{protID}"
		end
	end
	domain_tuples.each do |protID, doms|
		doms.each do |d|
			handler.puts "#{d}\t#{protID}"
		end
	end
	handler.close
end

def load_cath_info(cath_file, domain_class)
	proteins_2_domains = {}
	File.open(cath_file).each do |line|
		line.chomp!
		protID, superfam, funfam = line.split("\t")
		query = proteins_2_domains[protID]
		if domain_class == 'superfamilyID'
			domain_type = superfam
		else
			domain_type = funfam
		end
		if query.nil?
			proteins_2_domains[protID] = [domain_type]
		else
			query << domain_type
		end

	end
	return proteins_2_domains
end

def generate_tuples(prot_annots, prot_domains, go)
	domain_tuples = {}
	annot_tuples = {}
	term_levels = []
	unless go.nil?
		terms = prot_annots.values.flatten.uniq
		term_levels = terms.map{|a| go.get_term_level(a.to_sym)}
		min_level = term_levels.min
	end
	prot_annots.each do |protID, annots|
		domains = prot_domains[protID]
		unless domains.nil?
			domains.each do |d|	
				query = domain_tuples[protID]
				if query.nil?
					domain_tuples[protID] = [d]
				else
					query << d unless query.include?(d)
				end
			end
		end
		paths = []
		paths = annots.map{|a| go.get_parental_path(a.to_sym, :shortest_path, min_level)}.flatten if !go.nil?
		(annots | paths).each do |a|
			query = annot_tuples[protID]
			if query.nil?
				annot_tuples[protID] = [a]
			else
				query << a unless query.include?(a)
			end
		end
	end
	return domain_tuples, annot_tuples
end

##########################
#OPT-PARSER
##########################
options = {}
OptionParser.new do |opts|
	opts.banner = "Usage: #{__FILE__} [options]"

	options[:protein_domains] = nil
	opts.on("-a", "--protein_domains PATH", "Training proteins with CATH domains") do |data|
		options[:protein_domains] = data
	end

	options[:annotated_proteins] = nil
	opts.on("-b", "--annotated_proteins PATH", "Training proteins with annotations") do |data|
		options[:annotated_proteins] = data
	end

	options[:domain_class] = 'funfamID'
	opts.on("-d", "--domain_class STRING", "Domain identifiers type. Please choose funfamID or superfamilyID") do |data|
		options[:domain_class] = data
	end

	options[:gene_ontology] = nil
	opts.on("-g", "--gene_ontology PATH", "Gene ontology file") do |data|
		options[:gene_ontology] = data
	end

	options[:output_network] = 'tripartite_network.txt'
	opts.on("-o", "--output_network PATH", "Output tripartite network from CAFA information") do |data|
		options[:output_network] = data
	end

	opts.on_tail("-h", "--help", "Tool information") do
		puts opts
		exit
	end

end.parse!

##########################
#MAIN
##########################
go = nil
go = Ontology.new(file: options[:gene_ontology], load_file: true) unless options[:gene_ontology].nil?
proteins_with_annotations = load_hash(options[:annotated_proteins], 'a')
proteins_2_domains = load_cath_info(options[:protein_domains], options[:domain_class])
domain_tuples, annot_tuples = generate_tuples(proteins_with_annotations, proteins_2_domains, go)
build_tripartite_network(domain_tuples, annot_tuples, options[:output_network])