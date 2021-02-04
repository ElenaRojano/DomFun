#! /usr/bin/env ruby

##########################
# Rojano E. & Seoane P., Feb 2021
# Select curated annotations from UniProt
##########################

require 'optparse'

KEGG = 0
REACT = 1
GOMF = 2
GOBP = 3
GOCC = 4

##########################
#METHODS
##########################

def load_uniProt_file(input_file)
	uniProt_data = {}
	header = true
	File.open(input_file).each do |line|
		line.chomp!
		if !header
			uniProt_proteins_info = line.split("\t", 7)
			prot_id, fragment = uniProt_proteins_info.shift(2)
			next if fragment == 'fragment' || uniProt_proteins_info.all? {|a| a == ""}
			kegg, reactome, gomf, gobp, gocc = uniProt_proteins_info
			annots = [
				kegg.split(';'), 
				reactome.split(';'),
				gomf.scan(/GO:\d{7,7}/),
				gobp.scan(/GO:\d{7,7}/),
				gocc.scan(/GO:\d{7,7}/)
			]
			uniProt_data[prot_id] = annots
		else
			header = false
		end
	end
	return uniProt_data
end

def load_curated_file(file)
	curated_proteins_annots = {}
	File.open(file).each do |line|
		line.chomp!
		proteinID, annot = line.split("\t")
		query = curated_proteins_annots[proteinID]
		if query.nil?
			curated_proteins_annots[proteinID] = [annot]
		else
			query << annot
		end
	end
	return curated_proteins_annots
end

def select_curated_proteins(uniProt_data, curated_proteins_annots)
	curated_info = {}
	uniProt_data.each do |protID, annots|
		annots2save = []
		if annots[GOMF].empty? && annots[GOBP].empty? && annots[GOCC].empty?
			annots2save = annots
		else
			annots2save.concat(annots[KEGG..REACT])
			curated_annots = curated_proteins_annots[protID] #find curated GO
			if !curated_annots.nil? #protein with GO curated annotations
				[GOMF, GOBP, GOCC].each do |go|
					annots2save << annots[go] & curated_annots
				end
			else #if proteins have no GO curated, check if they have KEGG/Reactome
				annots2save.concat([[], [], []]) # GOMF, GOBP, GOCC
			end
		end
		curated_info[protID] = annots2save if !annots2save.all?{|a| a.empty?}
	end
	return curated_info
end


##########################
#OPT-PARSER
##########################
options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"

  options[:input_file] = nil
  opts.on("-a", "--input_file PATH", "Input file with protein annotations (from UniProt)") do |data|
    options[:input_file] = data
  end

  options[:curated_annots] = nil
  opts.on("-c", "--curated_annots PATH", "Input file with proteins and curated annotations") do |data|
    options[:curated_annots] = data
  end

  options[:output_file] = 'curated_proteins.txt'
  opts.on("-o", "--output_file PATH", "Output file with UniProt terms with GO terms curated") do |data|
    options[:output_file] = data
  end

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!

##########################
#MAIN
##########################
uniProt_data = load_uniProt_file(options[:input_file])
curated_proteins_annots = load_curated_file(options[:curated_annots])
final_curated_proteins = select_curated_proteins(uniProt_data, curated_proteins_annots)

File.open(options[:output_file], 'w') do |f|
	f.puts "Entry\tCross-reference (kegg)\tCross-reference (reactome)\tGene ontology (molecular function)\tGene ontology (biological process)\tGene ontology (cellular component)"
	final_curated_proteins.each do |protein, annotations|
		f.puts "#{protein}\t#{annotations.map{|a| a.join(';')}.join("\t")}"
	end
end