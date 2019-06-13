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
require 'generalMethods.rb'
require 'csv'
require 'optparse'
require "statistics2"
require "terminal-table"

# Domain to functional annotation predictor
# Based on domain-annotation association, this predictor can add functions to a group of domains of a protein
# It predict the most putative functions associated to a protein based on their domains.

##########################
#METHODS
##########################

def find_domains(cath_data, gene2proteins, protein_to_search)
  domain_ids = []
  protein_ids = gene2proteins[protein_to_search].uniq
  protein_ids.each do |protein_id|
    domain_ids << cath_data[protein_id]
  end
  return domain_ids.flatten.uniq
end




def load_domain_to_pathway_association(associations_file)
	domain_to_pathway_associations = {}
	File.open(associations_file).each do |line|
		line.chomp!
		annotation, domain, association_value = line.split("\t")
    association_value = association_value.to_f
		query = domain_to_pathway_associations[domain]
		if query.nil?
			domain_to_pathway_associations[domain] = [[annotation, association_value]]
		else
			query << [annotation, association_value]
		end
	end
	return domain_to_pathway_associations
end

def load_domains_to_predict(domains_file)
  domains_to_predict = {}
  File.open(domains_file).each do |line|
    line.chomp!
    protein_id, domains = line.split("\t")
    domains_to_predict[protein_id] = domains.split(',')
  end
  return domains_to_predict
end


      # #hyper must be ln not log10 from net analyzer
      # #https://en.wikipedia.org/wiki/Fisher%27s_method
      # lns = associations.map{|a| Math.log(10 ** -a)} #hyper values come as log10 values
      # sum = lns.inject(0){|s, a| s + a} 
      # combined_pvalue = Statistics2.chi2_x(sample_length *2, -2*sum)  
      # regionAttributes_array[i] << combined_pvalue


def search4function(domains_to_predict, domain_to_pathway_associations)
  domain_to_function_and_association_value = {}
  domains_to_predict.each do |domain|
    associations = domain_to_pathway_associations[domain]
    if !associations.nil?
      domain_to_function_and_association_value[domain] = associations
    end
  end
  return domain_to_function_and_association_value
end


def group_by_function(domain_to_function_and_association_value)
  function_to_domains = {}
  functionAttributes = {}
  association_scores = {}
  domain_to_function_and_association_value.each do |domain, annotations|
    annotations.each do |annotation_id, association_score|
      query = function_to_domains[annotation_id]
      if query.nil?
        function_to_domains[annotation_id] = [domain]
      else
        query << domain
      end
      functionAttributes[annotation_id] = []
      query = association_scores[annotation_id]
      if query.nil?
        association_scores[annotation_id] = {domain => association_score}
      else
        query[domain] = association_score 
      end
    end   
  end
  return function_to_domains, functionAttributes, association_scores
end

def generate_domain_annotation_matrix(function_to_domains, association_scores, domains_to_predict, null_value=0)
  # #method for creating the hpo to region matrix for plotting
  # #info2predict = hpo list from user
  # #hpo_associated_regions = [[chr, start, stop, [hpos_list], [weighted_association_scores]]]
  domain_annotation_matrix = []
  function_to_domains.each do |function_ID, domains_list|
    row = []
    domains_to_predict.each do |user_domain|
      value = association_scores[function_ID][user_domain]
      if value.nil?
        row << null_value
      else
        row << value
      end
    end
    domain_annotation_matrix << row
  end
  return domain_annotation_matrix
end

def scoring_regions(functionAttributes, domain_annotation_matrix, scoring_system='fisher', pvalue_cutoff=0, freedom_degree='maxnum', null_value=0)
  #hpo_associated_regions = [[chr, start, stop, [hpos_list], [weighted_association_scores]]]
  #hpo_region_matrix = [[0, 0.4, 0, 0.4], [0, 0, 0.5, 0.4]]
  functionAttributes_array = functionAttributes.values
  max_cluster_length = domain_annotation_matrix.map{|x| x.count {|i| i != 0}}.max if freedom_degree == 'maxnum'
  domain_annotation_matrix.each_with_index do |associations, i|
    sample_length = nil
    if freedom_degree == 'maxnum'
      sample_length = max_cluster_length
    else
      abort("Invalid freedom degree calculation method: #{freedom_degree}")
    end

    if scoring_system == 'fisher'
      #hyper must be ln not log10 from net analyzer
      #https://en.wikipedia.org/wiki/Fisher%27s_method
      lns = associations.map{|a| Math.log(10 ** -a)} #hyper values come as log10 values
      sum = lns.inject(0){|s, a| s + a} 
      combined_pvalue = Statistics2.chi2_x(sample_length *2, -2*sum)
      functionAttributes_array[i] << combined_pvalue
      # STDERR.puts functionAttributes.inspect
      # Process.exit 
    end
  end
end


##########################
#OPT-PARSER
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"
  
  options[:input_domains] = nil
  opts.on("-d", "--input_domains PATH", "Input file with protein domains to predict their function, separated by ','") do |data|
    options[:input_domains] = data
  end

  options[:protein_domains_file] = nil
  opts.on("-f", "--protein_domains_file PATH", "Input protein-domains file from CATH") do |data|
    options[:protein_domains_file] = data
  end

  options[:input_associations] = nil
  opts.on("-a", "--input_associations PATH", "Domain-function associations") do |data|
    options[:input_associations] = data
  end

  options[:domain_category] = "superfamilyID"
  opts.on("-c", "--domain_category PATH", "Domain category. Please choose one: superfamilyID or funfamID" ) do |data|
    options[:domain_category] = data
  end

  options[:input_protein] = "MACF1"
  opts.on("-p", "--input_protein PATH", "Protein to predict in UniProt ID" ) do |data|
    options[:input_protein] = data
  end


  # options[:best_thresold] = 0.5
  # opts.on("-b", "--best_thresold FLOAT", "Association value thresold") do |best_thresold|
  #   options[:best_thresold] = best_thresold.to_f
  # end

  options[:multiple_domains] = false
    opts.on("-u", "--multiple_domains", "Set if multiple profiles") do
  options[:multiple_domains] = true
  end

end.parse!

##########################
#MAIN
##########################

null_value = 0
cath_data, protein2gene, gene2proteins, cath_proteins_number = load_cath_data(options[:protein_domains_file], options[:domain_category])
domain_to_pathways_associations = load_domain_to_pathway_association(options[:input_associations])
domains_to_predict = find_domains(cath_data, gene2proteins, options[:input_protein])
Process.exit

domains_to_predict = load_domains_to_predict(options[:input_domains])
domains_to_predict.each do |protein, domains|
  domain_to_function_and_association_value = search4function(domains, domain_to_pathways_associations)
  function_to_domains, functionAttributes, association_scores = group_by_function(domain_to_function_and_association_value)
  annotation_matrix = generate_domain_annotation_matrix(function_to_domains, association_scores, domains, 0)
  scoring_regions(functionAttributes, annotation_matrix, 'fisher', 0, 'maxnum', null_value)
  #at this point, Fisher has been performed
  # STDERR.puts function_to_domains.inspect
  functionAttributes.each do |function, combined_pvalue|
    domains = function_to_domains[function]
    puts "#{domains.join(',')}\t#{function}\t#{combined_pvalue}"
  end
# association_scores => {"R-HSA-114608"=>{"1.10.287.600"=>0.6506696909373191, "1.10.246.10"=>2.9913972767050714, "1.10.238.10"=>0.6474149117916567, "1.10.150.300"=>1.879010877888625},...


end
  # domain_to_annotation, domainAttributes, association_scores = group_by_domain(domains_with_annotations)
  # annotation_matrix = generate_domain_annotation_matrix(domain_to_annotation, association_scores, domains, 0)
  # scoring_regions(domainAttributes, annotation_matrix, 'fisher', 0, 'maxnum', null_value)
  # association_scores.each do |annotation, associations|
  #   #printar salida de la tabla, comprobar por que se repite muchas veces (mal metido en el hash)
  #   STDERR.puts association_scores.inspect
  #   # STDERR.puts association_values
  # end


#retrieve_pathways_from_domains(domain_to_pathways_associations, domains_to_predict)
