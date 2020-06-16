#! /usr/bin/env ruby
##########################
# Rojano E. & Seoane P., July 2019
# Normalize combined scores for their use with CAFA validation
##########################

require 'optparse'

##########################
#METHODS
##########################
module Enumerable

    def sum
      self.inject(0){|accum, i| accum + i }
    end

    def mean
      self.sum/self.length.to_f
    end

    def sample_variance
      m = self.mean
      sum = self.inject(0){|accum, i| accum +(i-m)**2 }
      sum/(self.length - 1).to_f
    end

    def standard_deviation
      Math.sqrt(self.sample_variance)
    end

end

def load_predictions(input_file)
	predictions_data = []
	File.open(input_file).each do |line|
		line.chomp!
		protID, domains, funSys, combScore = line.split("\t")
		predictions_data << [protID, domains, funSys, combScore.to_f]
	end
	return predictions_data
end

##########################
#OPT-PARSER
##########################
options = {}
OptionParser.new do |opts|
	opts.banner = "Usage: #{__FILE__} [options]"

	options[:input_file] = nil
		opts.on("-a", "--input_file PATH", "Input file with association values to normalize") do |data|
		options[:input_file] = data
	end

	options[:integration_method] = 'fisher'
		opts.on("-i", "--integration_method STRING", "Integration method") do |data|
		options[:integration_method] = data
	end

	options[:normalization_mode] = 'normal'
		opts.on("-m", "--normalization_mode STRING", "Normalization mode: normal or max") do |data|
		options[:normalization_mode] = data
	end

	options[:output_file] = 'normalized_associations.txt'
		opts.on("-o", "--output_file PATH", "Output association file with normalized values") do |data|
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
predictions_data = load_predictions(options[:input_file])
standardized_data = []

if options[:integration_method] == 'fisher'
	predictions_data.each do |protID, domains, funSys, combScore|
		stdScore = 1 - combScore
		#CAFA validation score must be in range (0,1]
		#Remember that fisher's integration most significant value is 0
		standardized_data << [protID, domains, funSys, stdScore] if stdScore > 0.001  
	end
else
	#https://www.researchgate.net/post/How_do_i_normalize_data_from_0_to_1_range
	combScores = predictions_data.map{|a| a[3] }
	combScoresAverage = combScores.mean
	combScoresSD = combScores.standard_deviation
	maxCombScore = combScores.max
	predictions_data.each do |protID, domains, funSys, combScore|
		if options[:normalization_mode] == 'normal'
			score = (combScore - combScoresAverage).fdiv(combScoresSD)
			if score > 2
				score = 2
			elsif score < -2
				score = -2
			end
			stdScore = score.fdiv(4) + 0.5	
		elsif options[:normalization_mode] == 'max'			
			stdScore = combScore.fdiv(maxCombScore)
		end				
		#CAFA validation score must be in range (0,1]
		standardized_data << [protID, domains, funSys, stdScore] #if stdScore > 0.001  
	end
end

handler = File.open(options[:output_file], 'w')
standardized_data.each do |data|
	handler.puts data.join("\t")
end
handler.close