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

##########################
#METHODS
##########################
def build_tripartite_networks(nomenclature_annotations, cath_data, path, protein2gene, translate2gene)
	records = Hash.new(0)
	#STDERR.puts cath_data.inspect
  nomenclature_annotations.each do |nomenclature, protein_annotations|
		annots = []
		datas = []
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
cath_data, protein2gene, cath_proteins_number, cath_data_supp = load_cath_data(options[:input_domains], options[:category_type])
nomenclature_annotations, number_of_proteins = load_proteins_file(options[:input_annotations], options[:annotation_types])

networks_path = nil
if options[:category_type] == 'funfamID'
	networks_path = 'PPP_results/networks/funfam_networks'
else
	networks_path = 'PPP_results/networks/superfamily_networks'
end
FileUtils.mkdir_p networks_path
puts "Generating tripartite networks. This can take a while, please wait."
protein_stats = build_tripartite_networks(nomenclature_annotations, cath_data, networks_path, protein2gene, options[:translate2gene])
File.open(options[:output_stats], 'w') do |f|
  protein_stats.each do |annotation_type, number_of_proteins|
  	f.puts "#{annotation_type}\t#{number_of_proteins}"
  end
  f.puts "Total of Uniprot proteins\t#{number_of_proteins}"
  f.puts "Total of CATH proteins\t#{cath_proteins_number}"
end
