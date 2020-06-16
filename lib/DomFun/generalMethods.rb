def load_proteins_file(file, annotation_types)
	protein_annotations = {}
	proteins_without_annotations = []
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
		fields = line.split("\t", fields_to_split + 1)
		protID = fields.shift
		annotation_types.each_with_index do |type, i|
			annotations = fields[i].split(/[;,]/)
			if !annotations.empty?
				if type.include?('go')
					go_annotations = []
					annotations.each do |go_term|
						go_name, go_id = go_term.split('GO:')
						go_annotations << "GO:".concat(go_id.tr(']', '')) unless go_id.nil?	
					end
					protein_annotations[type][protID] = go_annotations
				else
					protein_annotations[type][protID] = annotations
				end
			end
			if fields.count("") == fields_to_split
				proteins_without_annotations << protID
			end
		end
		counter += 1
	end
	return protein_annotations, counter, proteins_without_annotations.uniq
end

def load_cath_data(file, category, dictionary_key='gene_name')
	if dictionary_key == 'gene_name'
		field = 3
	elsif dictionary_key == 'geneID' # UNIPROT entry_name
		field = 4
	end
	cath_data = {}
	protein2gene_dict = {}
	csv_file = CSV.read(file, { :col_sep => "\t" })
	csv_file.delete_at(0)
	csv_file.each do |protein_domains_data|
		next if protein_domains_data.empty?
		protein_id = protein_domains_data[0]
		protein_alternative_name = protein_domains_data[field]
		next if protein_domains_data[3].include?('fusion') # Only can checked in cath gene name field
		protein_alternative_name.gsub!(' ', '_') if protein_alternative_name.include?(' ')
		superfamilyID = protein_domains_data[5]
		funfamID = protein_domains_data[6]
		term2save = nil
		if category == 'superfamilyID'
			term2save = superfamilyID
		elsif category == 'funfamID'
			term2save = funfamID
		end
		add_term2dictionary(cath_data, protein_id, term2save)
		protein2gene_dict[protein_id] = protein_alternative_name if protein_alternative_name != 'NULL'
	end
	cath_proteins_number = cath_data.keys.length
	return cath_data, protein2gene_dict, cath_proteins_number
end

def invert_hash(hash)
	new_hash = {}
	hash.each do |k, v|
		query = new_hash[v]
		if query.nil?
			new_hash[v] = [k]
		else
			query << k
		end
	end
	return new_hash
end

def add_term2dictionary(dict, key, term)
	query = dict[key]
	if query.nil?
        	dict[key] = [term]
	else
		query << term
	end
end

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