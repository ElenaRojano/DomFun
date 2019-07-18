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
require 'generalMethods.rb'
require 'csv'
require 'optparse'
require "statistics2"
require "terminal-table"
require 'report_html'


##########################
#METHODS
##########################

def get_protein_domains(cath_data, protein)
  # puts cath_data.inspect
  # Process.exit
  unless cath_data[protein].nil?
    domains = cath_data[protein].uniq
    domains_to_predict = [protein, domains] #unless cath_data[protein].nil?
  end
  # STDERR.puts protein
  # STDERR.puts domains_to_predict
  return domains_to_predict
end


def load_domain_to_pathway_association(associations_file, threshold)
	domain_to_pathway_associations = {}
	File.open(associations_file).each do |line|
		line.chomp!
		annotation, domain, association_value = line.split("\t")
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

def scoring_regions(function_to_domains, domain_annotation_matrix, scoring_system='fisher', freedom_degree='maxnum', null_value=0)
  function_to_domains_array = function_to_domains.values
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
      function_to_domains_array[i] << combined_pvalue
    end
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

  options[:limit_show] = 0
  opts.on("-l", "--limit_show INTEGER", "Maximum number of predictions to show") do |limit_show|
    options[:limit_show] = limit_show.to_i
  end

  options[:output_file] = 'predictions_file.txt'
  opts.on("-o", "--output_file PATH", "Predictions file") do |data|
    options[:output_file] = data
  end

  options[:proteins_2predict] = nil
  opts.on("-p", "--proteins_2predict PATH", "Protein to predict. Please use UniProt IDs" ) do |data|
    options[:proteins_2predict] = data
  end

  options[:pvalue_threshold] = 0.05
  opts.on("-t", "--pvalue_threshold FLOAT", "P-value threshold") do |pvalue_threshold|
    options[:pvalue_threshold] = best_threshold.to_f
  end

  options[:association_threshold] = 2
  opts.on("-T", "--association_threshold FLOAT", "Association value threshold") do |association_threshold|
    options[:association_threshold] = association_threshold.to_f
  end


  options[:html_file] = "DomFunPredictions.html"
  opts.on("-R", "--html_file PATH", "HTML file with prediction summary") do |html_file|
    options[:html_file] = html_file
  end

  options[:html_reporting] = false
    opts.on("-r", "--html_reporting", "Set to generate an HTML report") do
  options[:html_reporting] = true
  end

  options[:unknown_proteins] = 'unknown_proteins_list.txt'
  opts.on("-U", "--unknown_proteins PATH", "List of proteins not found in the system") do |data|
    options[:unknown_proteins] = data
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
# 1. Load associations and domains data
cath_data, protein2gene, gene2proteins, cath_proteins_number = load_cath_data(options[:protein_domains_file], options[:domain_category])
# 2. Load protein(s) to predict
if File.exist?(options[:proteins_2predict])
  if !options[:multiple_proteins]
    options[:proteins_2predict] = [File.open(options[:proteins_2predict]).readlines.map!{|line| line.chomp}]
  else
    multiple_profiles = []
    File.open(options[:proteins_2predict]).each do |line|
      line.chomp!
      multiple_profiles << line.split('|')
    end
    options[:proteins_2predict] = multiple_profiles
  end
else
  if !options[:multiple_proteins]
    options[:proteins_2predict] = [options[:proteins_2predict].split('|')]
  else
    options[:proteins_2predict] = options[:proteins_2predict].split('!').map{|profile| profile.split('|')}
  end
end


# 3. Load domain-FunSys associations
domain_to_pathways_associations = load_domain_to_pathway_association(options[:input_associations], options[:association_threshold])

# 4. Prediction performance
domains_per_protein = {}
domains_without_FunSys = {}
options[:proteins_2predict].each do |protein|
  # STDERR.puts protein.join('').inspect
  domains_to_predict = get_protein_domains(cath_data, protein.join(''))
  unless domains_to_predict.nil?
    null_value = 0
    protein_id = domains_to_predict.shift
    domains_to_predict.each do |domains|
      # STDERR.puts domains.inspect
      # STDERR.puts data.inspect


      # protein_id = domains.shift
      # if domains.empty?
      #   handler = File.open(options[:unknown_proteins], 'w')
      #   handler.puts "Protein #{protein_id} has no domains found"
      #   handler.close
      # else     
        domain_function_assocValue = search4function(domains, domain_to_pathways_associations)
        # STDERR.puts domain_function_assocValue.inspect
        domains_without_FunSys[protein_id] = domains if domain_function_assocValue.empty?
          # handler = File.open(options[:output_file], 'w')
          # handler.puts "Protein #{protein_id} with domains #{domain(s).join(', ')} were not found to be associated with any FunSys"
          # handler.close
        
          #STDERR.puts domain_function_assocValue.inspect
        function_to_domains, association_scores = group_by_function(domain_function_assocValue)
        annotation_matrix = generate_domain_annotation_matrix(function_to_domains, association_scores, domains, 0)
        scoring_regions(function_to_domains, annotation_matrix, 'fisher', 'maxnum', null_value)
        predictions = []
        # STDERR.puts function_to_domains.inspect
        domains_group = []
        domains_association_values = []
        function_to_domains.each do |function, info|
          combined_score = info.pop
          function_domains_associated_and_value = association_scores[function]
          info.each do |domain|
            domains_association_values << function_domains_associated_and_value[domain]
            domains_group << domain
            # STDERR.puts domains_association_values.inspect
          end
          predictions << [protein_id, domains_group.join(', '), function, domains_association_values.map{|a| a.round(3)}.join(', '), combined_score] 
          # STDERR.puts predictions.inspect
          domains_association_values = []
          domains_group = []
        end
        predictions.sort!{|r1, r2| r1.last <=> r2.last}
        predictions.each do |info|
          STDERR.puts info.inspect
          puts info.join("\t")
        end
        # STDERR.puts predictions.inspect
        # handler = File.open(options[:output_file], 'w')
        # handler.puts "ProteinID\tDomains\tFunSys term\tAssociation values\tCombined score (P-value)"
        # if options[:limit_show] == 0
        #   predictions.each do |info|
        #     handler.puts info.join("\t") if info.last <= options[:pvalue_thresold]
        #   end
        # else
        #   predictions.each_with_index do |info, limit|
        #     handler.puts info.join("\t") if info.last <= options[:pvalue_thresold]
        #     break if limit == options[:limit_show] - 1
        #   end
        # end
        # domains_without_FunSys.each do |protein, domains|
        #   handler.puts "#{protein}\t#{domains.join(', ')}\t#{'-'}\t#{'-'}\t#{'-'}"
        # end
        # handler.close
        # STDERR.puts predictions.inspect
        report_data(predictions, options[:html_file]) if options[:html_reporting]
      # end
    end
  end
end



