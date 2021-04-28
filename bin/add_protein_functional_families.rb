#! /usr/bin/env ruby
##########################
# Rojano E. & Seoane P., June 2019
# Generate tripartite networks with domains-proteins-FunSys data
# Protein IDs and FunSys (GO-MF, KEGG and Reactome) from UniProtKB.
# Protein domains (Superfamilies and FunFams) from CATH.
##########################

REPORT_FOLDER=File.expand_path(File.join(File.dirname(__FILE__), '..', 'templates'))
ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'DomFun'))
require 'generalMethods'
require 'csv'
require 'optparse'
require 'fileutils'
require 'semtools'

##########################
#METHODS
##########################
def build_tripartite_networks(nomenclature_annotations, cath_data, path, protein2gene, translate2gene, go)
	records = Hash.new(0)
  nomenclature_annotations.each do |nomenclature, protein_annotations|
		annots = []
		datas = []
    if nomenclature == 'go' && !ontology.nil?
      terms = protein_annotations.values.flatten.uniq
      term_levels = terms.map{|a| go.get_term_level(a.to_sym)}
      min_level = term_levels.min
    end
		protein_annotations.each do |protID, annotations|
			query_cath_data = cath_data[protID]
      if !query_cath_data.nil?
        if !translate2gene
          recordID = protID
        else
  				recordID = protein2gene[protID]
  				recordID = protID if recordID.nil?
        end
				annotations.each do |annotation|
          if nomenclature == 'reactome'
            annotation.gsub!(/R-([A-Z]{3,3})-/, 'R-')
          end
					annots << [annotation, recordID]
				end
				query_cath_data.each do |data|
					datas << [data, recordID]
				end
			end
		end
		records[nomenclature] += annots.map{|pair| pair.last}.uniq.length
		File.open(File.join(path, "network_#{nomenclature}"), 'w') do |f|
			annots.uniq.each do |pair|
				f.puts pair.join("\t")
			end
			datas.uniq.each do |pair|
				f.puts pair.join("\t")
			end
		end
	end
	return records
end

def load_proteins_file(file, annotation_types)
  protein_annotations = {}
  annotation_types.each do |type| # initialize annotation hashes
    protein_annotations[type] = {}
  end
  fields_to_split = annotation_types.length 
  counter = 0
  File.open(file).each do |line|
    line.chomp!
    if counter == 0
      counter += 1
      next
    end
    line.gsub!(' ', '')
    fields = line.split("\t", fields_to_split + 2) #Sum according to the total of GO subontologies -1
    protID = fields.shift.to_sym
    annotation_types.each_with_index do |type, i|
      annotations = fields[i].split(/[;]/)
      if !annotations.empty?
        if type.include?('go')
          go_annotations = []
          annotations.each do |go_term|
            go_annotations << go_term unless go_term.nil? 
          end
          protein_annotations[type][protID] = go_annotations
        else
          protein_annotations[type][protID] = annotations
        end
      end
    end
    counter += 1
  end
  return protein_annotations, counter
end


##########################
#OPT-PARSER
##########################
options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:input_annotations] = nil
  opts.on("-a", "--input_annotations PATH", "Input file with gene annotations") do |data|
    options[:input_annotations] = data
  end

  options[:calculate_proteins_by_domain] = false
  opts.on("-c", "--calculate_proteins_by_domain", "Calculate the number of proteins that a domain has") do
    options[:calculate_proteins_by_domain] = true
  end

  options[:input_domains] = nil
  opts.on("-d", "--input_domains PATH", "Input file with protein domains") do |data|
    options[:input_domains] = data
  end

  options[:search_domain] = true
  opts.on("-f", "--search_domain", "Search full protein domains. If false, search funfams") do
    options[:search_domain] = false
  end

  options[:gene_ontology] = nil
  opts.on("-g", "--gene_ontology PATH", "Gene ontology file") do |data|
    options[:gene_ontology] = data
  end

  options[:save_network] = nil
  opts.on("-n", "--save_network PATH", "Path to save network files") do |data|
    options[:save_network] = data
  end

  options[:output_file] = 'uniprot_translated.txt'
  opts.on("-o", "--output_file PATH", "Output file with UniProt GeneName structure for prediction") do |data|
    options[:output_file] = data
  end

  options[:annotation_types] = %w[ kegg reactome go]
  opts.on("-p", "--annotation_types STRING", "List of annotation types separated by commas") do |data|
    options[:annotation_types] = data.split(",")
  end

  options[:output_stats] = 'uniprot_stats.txt'
  opts.on("-s", "--output_stats PATH", "Output file with UniProt stats") do |data|
    options[:output_stats] = data
  end

  options[:translate2gene] = false
  opts.on("-T", "--translate2gene", "Translate proteins to genes") do
    options[:translate2gene] = true
  end

  options[:category_type] = 'funfamID'
  opts.on("-t", "--category_type STRING", "Input category of domains. Options: funfamID, superfamilyID") do |data|
    options[:category_type] = data
  end  

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!

##########################
#MAIN
##########################

puts "Loading data..."
go = nil
go = Ontology.new(file: options[:gene_ontology], load_file: true) unless options[:gene_ontology].nil?
cath_data, protein2gene, cath_proteins_number, cath_data_supp = load_cath_data(options[:input_domains], options[:category_type])
nomenclature_annotations, number_of_proteins = load_proteins_file(options[:input_annotations], options[:annotation_types])

networks_path = nil
if options[:category_type] == 'funfamID'
	networks_path = File.join(options[:save_network], 'funfam_networks')
else
	networks_path = File.join(options[:save_network], 'superfamily_networks')
end
FileUtils.mkdir_p networks_path
puts "Generating tripartite networks. This can take a while, please wait."
protein_stats = build_tripartite_networks(nomenclature_annotations, cath_data, networks_path, protein2gene, options[:translate2gene], go)
File.open(options[:output_stats], 'w') do |f|
  protein_stats.each do |annotation_type, number_of_proteins|
  	f.puts "#{annotation_type}\t#{number_of_proteins}"
  end
  f.puts "Total of Uniprot proteins\t#{number_of_proteins}"
  f.puts "Total of CATH proteins\t#{cath_proteins_number}"
end
