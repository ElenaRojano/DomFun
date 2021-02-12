#! /usr/bin/env ruby
##########################
# Rojano E. & Seoane P., June 2019
# Domain to functional annotation predictor
# Based on domain-annotation association, this predictor can add functions to a group of domains of a protein
# It predict the most putative functions associated to a protein based on their domains.
# Protein IDs and FunSys (GO-MF, KEGG and Reactome) from UniProtKB.
# Protein domains (Superfamilies and FunFams) from CATH.
##########################

REPORT_FOLDER=File.expand_path(File.join(File.dirname(__FILE__), '..', 'templates'))
ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'DomFun'))
require 'generalMethods'
require 'csv'
require 'optparse'
require "statistics2"
require "terminal-table"
require 'report_html'
require 'parallel'

##########################
#METHODS
##########################

def get_protein_domains(cath_data, protein, gene2protein, identifier_mode)
  domains_to_predict = []
  if identifier_mode == 'mixed'
    proteins = gene2protein[protein]
    if !proteins.nil?
      proteins.each do |protein|
        domains_to_predict.concat(cath_data[protein])
      end
    else
      domains_to_predict = cath_data[protein]
    end
  else
    domains_to_predict = cath_data[protein]
  end
  domains_to_predict = [] if domains_to_predict.nil?
  return domains_to_predict.uniq
end

def load_domain_to_pathway_association(associations_file, threshold, white_list=nil)
	domain_to_pathway_associations = {}
	File.open(associations_file).each do |line|
		line.chomp!
		annotation, domain, association_value = line.split("\t")
    next if !white_list.nil? && white_list[domain].nil?
    association_value = association_value.to_f
		next if association_value < threshold
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
  association_scores = {}
  domain_to_function_and_association_value.each do |domain, annotations|
    annotations.each do |annotation_id, association_score|
      query = function_to_domains[annotation_id]
      if query.nil?
        function_to_domains[annotation_id] = [domain]
      else
        query << domain
      end
      query = association_scores[annotation_id]
      if query.nil?
        association_scores[annotation_id] = {domain => association_score}
      else
        query[domain] = association_score 
      end
    end   
  end
  return function_to_domains, association_scores
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

def scoring_funsys(function_to_domains, domain_annotation_matrix, scoring_system, freedom_degree='maxnum', null_value=0, pvalue_threshold)
  domains_array = function_to_domains.values
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
      domains_array[i] << combined_pvalue
    elsif scoring_system == 'harmonic'
      lns = associations.map{|a| 10 ** -a}
      inv = lns.map{|n| 1.fdiv(n)}
      sum = inv.inject(0){|s,x| s + x}
      combined_pvalue = associations.length.fdiv(sum)
      domains_array[i] << combined_pvalue
    elsif scoring_system == 'stouffer'
      sum = associations.inject(0){|s,x| s + x}      
      combined_z_score = sum/Math.sqrt(sample_length)
      domains_array[i] << combined_z_score
    elsif scoring_system == 'average'
      sum = associations.inject(0){|s,x| s + x.abs}.fdiv(associations.length)
      domains_array[i] << sum
    elsif scoring_system == 'sum'
      sum = associations.inject(0){|s,x| s + x.abs}
      domains_array[i] << sum
    else
      abort("Invalid integration method: #{scoring_system}")
    end
  end
  if scoring_system == 'fisher' || scoring_system == 'harmonic'
    function_to_domains.select!{|function, attributes| attributes.last <= pvalue_threshold}
  else
    function_to_domains.select!{|function, attributes| attributes.last >= pvalue_threshold}
  end
end


def report_data(predictions, html_file)
  container = {:predictions => predictions }
  template = File.open(File.join(REPORT_FOLDER, 'report_data.erb')).read
  report = Report_html.new(container, 'Protein domains and FunSys predictions summary')
  report.build(template)
  report.write(html_file)
end

##########################
#OPT-PARSER
##########################

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"
  
  options[:input_associations] = nil
  opts.on("-a", "--input_associations PATH", "Domain-function associations") do |data|
    options[:input_associations] = data
  end

  options[:domain_category] = "superfamilyID"
  opts.on("-c", "--domain_category PATH", "Domain category. Please choose one: superfamilyID or funfamID" ) do |data|
    options[:domain_category] = data
  end

  options[:protein_domains_file] = nil
  opts.on("-f", "--protein_domains_file PATH", "Input protein-domains file from CATH") do |data|
    options[:protein_domains_file] = data
  end

  options[:integration_method] = 'fisher'
  opts.on("-i", "--integration_method STRING", "Integration method") do |data|
    options[:integration_method] = data
  end

  options[:identifier_mode] = 'normal'
  opts.on("-I", "--identifier_mode STRING", "Identifier mode: normal or mixed") do |data|
    options[:identifier_mode] = data
  end

  options[:output_file] = 'predictions_file.txt'
  opts.on("-o", "--output_file PATH", "Predictions file") do |data|
    options[:output_file] = data
  end

  options[:proteins_2predict] = nil
  opts.on("-p", "--proteins_2predict PATH", "Protein to predict. Please use UniProt IDs" ) do |data|
    options[:proteins_2predict] = data
  end

  options[:threads] = 0
  opts.on("-P", "--threads INTEGER", "Number of threads to parallelize") do |data|
    options[:threads] = data.to_i - 1
  end

  options[:pvalue_threshold] = 0.05
  opts.on("-t", "--pvalue_threshold FLOAT", "P-value threshold") do |pvalue_threshold|
    options[:pvalue_threshold] = pvalue_threshold.to_f
  end

  options[:association_threshold] = 0
  opts.on("-T", "--association_threshold FLOAT", "Association value threshold") do |association_threshold|
    options[:association_threshold] = association_threshold.to_f
  end

  options[:multiple_proteins] = false
    opts.on("-u", "--multiple_proteins", "Set if multiple profiles") do
  options[:multiple_proteins] = true
  end

  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end

end.parse!

##########################
#MAIN
##########################

# 1. Load protein(s) to predict
if File.exist?(options[:proteins_2predict])
  if !options[:multiple_proteins]
    options[:proteins_2predict] = [File.open(options[:proteins_2predict]).readlines.map!{|line| line.chomp.to_sym}]
  else
    multiple_proteins = []
    File.open(options[:proteins_2predict]).each do |line|
      multiple_proteins << line.chomp.to_sym
    end
    options[:proteins_2predict] = multiple_proteins
  end
else
  if !options[:multiple_proteins]
    options[:proteins_2predict] = [options[:proteins_2predict].split('|').map{|pt| pt.to_sym}]
  else
    options[:proteins_2predict] = options[:proteins_2predict].split('!').map{|profile| profile.split('|').map{|pt| pt.to_sym}}
  end
end

# 2. Load protein domains classification to get domains from proteins to predict
pt_white_list = {}
options[:proteins_2predict].each do |pt|
  pt_white_list[pt] = true
end
cath_data, protein2gene, cath_proteins_number = load_cath_data(options[:protein_domains_file], options[:domain_category], pt_white_list)
pt_white_list = nil

# 3. Load domain-FunSys associations
dm_white_list = {}
cath_data.each do |pt, domains|
  domains.each do |dm|
    dm_white_list[dm] = true
  end
end
domain_to_pathways_associations = load_domain_to_pathway_association(options[:input_associations], options[:association_threshold], dm_white_list)
dm_white_list = nil

# 4. Prediction
#handler = File.open(options[:output_file], 'w')
gene2protein = invert_hash(protein2gene) if options[:identifier_mode] == 'mixed'
all_predictions = Parallel.map(options[:proteins_2predict], in_process: options[:threads]) do |protein|
  domains = get_protein_domains(cath_data, protein, gene2protein, options[:identifier_mode])
  if domains.empty?
    nil
  else
    null_value = 0
    domain_function_assocValue = search4function(domains, domain_to_pathways_associations)
    function_to_domains, association_scores = group_by_function(domain_function_assocValue) 
    annotation_matrix = generate_domain_annotation_matrix(function_to_domains, association_scores, domains, 0) 
    
    scoring_funsys(
      function_to_domains, 
      annotation_matrix, 
      options[:integration_method], 
      'maxnum', 
      null_value, 
      options[:pvalue_threshold]
      )

      final_predictions = []
      function_to_domains.each do |funsys, domains_data|
        score = domains_data.pop
        #handler.puts "#{protein}\t#{domains_data.join(',')}\t#{funsys}\t#{score}"
        final_predictions << [protein.to_s, domains_data.join(','), funsys, score]
      end
      final_predictions
    end
end
all_predictions.each do |protein_predictions|
  unless protein_predictions.nil?
    protein_predictions.each do |info|
      puts info.join("\t")
    end
  end
end
 #<< "#{protein}\t#{domains_data.join(',')}\t#{funsys}\t#{score}"


#handler.close