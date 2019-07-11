#! /usr/bin/env ruby
##########################
# Rojano E. & Seoane P., July 2019
# Generate GO tripartite networks filtered with CAFA data
# For its use as system control (only GO)
##########################

REPORT_FOLDER=File.expand_path(File.join(File.dirname(__FILE__), '..', 'templates'))
ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'DomFun'))
require 'optparse'

##########################
#METHODS
##########################

def load_cafa_data(cafa_file)
	cafa_data = {}
	File.open(cafa_file).each do |line|
		line.chomp!
		next if line.include?('GO_Ont')
		cafa_info = line.split("\t")
		next unless cafa_info[1] == 'MF'
		go_term = cafa_info[4]
		gene_name = cafa_info[6]
		next if gene_name == 'NA'
		query = cafa_data[gene_name]
		if query.nil?
			cafa_data[gene_name] = [go_term]
		else
			query << go_term
		end
	end
	return cafa_data
end

def load_network_data(network)
	go_gene_rels = {}
	domain_gene_rels = {}
	File.open(network).each do |line|
		line.chomp!
		term, gene = line.split("\t")
		if term.include?('GO:')
			query = go_gene_rels[gene]
			if query.nil?
				go_gene_rels[gene] = [term]
			else
				query << term
			end
		else
			query = domain_gene_rels[gene]
			if query.nil?
				domain_gene_rels[gene] = [term]
			else
				query << term
			end
		end
	end
	return go_gene_rels, domain_gene_rels
end

def check_genes(term_protein_rels, cafa_data)
	cafa_data.each do |gene, go_list|
		term_protein_rels.delete(gene)
	end
	return term_protein_rels
end
##########################
#OPT-PARSER
##########################
options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:cafa_data] = nil
  opts.on("-a", "--cafa_data PATH", "Input CAFA gene annotations") do |data|
    options[:cafa_data] = data
  end

  options[:input_network] = nil
  opts.on("-n", "--input_network PATH", "Input network to parse") do |data|
    options[:input_network] = data
  end

  options[:output_network] = 'output_network.txt'
  opts.on("-o", "--output_network PATH", "Input network to parse") do |data|
    options[:output_network] = data
  end

  opts.on_tail("-h", "--help", "Show tool help") do
  	puts opts
  	exit
  end

end.parse!

##########################
#MAIN
##########################
cafa_data = load_cafa_data(options[:cafa_data])
go_gene_rels, domain_gene_rels = load_network_data(options[:input_network])
genes2gos_layer = check_genes(go_gene_rels, cafa_data)
genes2domains_layer = check_genes(domain_gene_rels, cafa_data)
handler = File.open(options[:output_network], 'w')
genes2gos_layer.each do |gene, go_list|
	go_list.each do |go_term|
		handler.puts "#{go_term}\t#{gene}"
	end
end
genes2domains_layer.each do |gene, domains|
	domains.each do |domain|
		handler.puts "#{domain}\t#{gene}"
	end
end
handler.close