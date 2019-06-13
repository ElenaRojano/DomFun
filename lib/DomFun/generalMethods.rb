def load_cath_data(file, category)
	cath_data = {}
	protein2gene = {}
	gene2proteins = {}
	csv_file = CSV.read(file, { :col_sep => "\t" })
	csv_file.delete_at(0)
	csv_file.each do |protein_domains_data|
		next if protein_domains_data.empty?
		protein_id = protein_domains_data[0]
		gene_name = protein_domains_data[3]
		next if gene_name.include?('fusion')
		gene_name = gene_name.gsub(' ', '_') if gene_name.include?(' ')
		superfamilyID = protein_domains_data[5]
		funfamID = protein_domains_data[6]
		term2save = nil
		if category == 'superfamilyID'
			term2save = superfamilyID
		elsif category == 'funfamID'
			term2save = funfamID
		end
		add_term2dictionary(cath_data, protein_id, term2save)
		protein2gene[protein_id] = gene_name if gene_name != 'NULL'
		query = gene2proteins[gene_name]
	    if query.nil?
	      	gene2proteins[gene_name] = [protein_id] if protein_id != 'NULL'
	    else
	      query << protein_id if protein_id != 'NULL'
		end
	end
	cath_proteins_number = cath_data.keys.length
	return cath_data, protein2gene, gene2proteins, cath_proteins_number
end

def add_term2dictionary(dict, key, term)
	query = dict[key]
	if query.nil?
        	dict[key] = [term]
	else
		query << term
	end
end
