#! /usr/bin/env ruby
##########################
# Rojano E. & Seoane P., July 2019
# Generate control from CAFA data (proteins and GO terms)
##########################

REPORT_FOLDER=File.expand_path(File.join(File.dirname(__FILE__), '..', 'templates'))
ROOT_PATH = File.dirname(__FILE__)
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'DomFun'))
require 'optparse'
require 'generalMethods.rb'

##########################
#OPT-PARSER
##########################
options = {}
OptionParser.new do |opts|
	opts.banner = "Usage: #{__FILE__} [options]"

	options[:input_cafa] = nil
		opts.on("-a", "--input_cafa PATH", "Input CAFA file") do |data|
		options[:input_cafa] = data
	end

	options[:output_file] = 'cafa_control.txt'
		opts.on("-o", "--output_file PATH", "Output control from CAFA data") do |data|
		options[:output_file] = data
	end

	opts.on_tail("-h", "--help", "Tool information") do
		puts opts
		exit
	end

end.parse!

##########################
#MAIN
##########################
cafa_data = load_cafa_data(options[:input_cafa])
handler = File.open(options[:output_file], 'w')
cafa_data.each do |protein, goTerms|
	handler.puts "#{protein}\t#{goTerms.join(';')}"
end
handler.close