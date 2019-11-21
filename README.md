# DomFun

DomFun is a new system to assign functions to unknown proteins using a systemic approach without considering their sequence but their domains associated with functional systems. It uses associations calculated between protein domains and functional annotations as training dataset and performs predictions over proteins (using UniProt identifiers) by finding their domains and if they have been associated with functional annotations (in GO molecular functions, biological processes, KEGG and Reactome pathway terms). 

## Installation

Add this line to your application's Gemfile:

```ruby
gem 'DomFun'
```

And then execute:

    $ bundle

Or install it yourself as:

    $ gem install DomFun

## Usage

To use DomFun, it is necessary to calculate associations between protein domains and functional annotations. For this, use NetAnalyzer ruby gem (https://rubygems.org/gems/NetAnalyzer), choosing the association index that fits the best to your data. Once these associations are calculated, they will be used to train DomFun and predict a set of proteins (giving UniProt identifiers) with unknown function.

Procedure steps performed in the study to develop DomFun follows:

    1. Generate relationships between domains, proteins and FunSys.
    2. Generate networks and calculate associations between domains and FunSys. 
    3. Training DomFun with domain-FunSys associations.

### STEP 1. Generate  domains-proteins-FunSys relationships ###
To generate the domains-FunSys associations to train DomFun is necessary to have, on the one hand proteins (UniProt identifiers preferably) with functional annotations (GO, KEGG or Reactome, for example) and, on the other hand, proteins with their functional domains classified (CATH or SCOP, for example). If domains-proteins-FunSys have been previously established this procedure will not be necessary.

Before calculating the association between protein domains and functional annotations/systems (abbr. FunSys), is necessary to construct the relationships between domains, proteins and FunSys:
*Note: DomFun was developed using the human proteins dataset; however it can be used in any species if protein identifiers, their domains and functional annotations are available.*
```
    add_protein_functional_families.rb -a proteinAnnotationsFile -d cathDomainsFile -p activeAnnots
``` 
Where
* -a proteinAnnotationsFile = Protein annotations file 
* -d cathDomainsFile = Domains file
* -p activeAnnots = List with annotations names used (KEGG, GO, Reactome...)

### STEP 2. Create tripartite network and association calculation ###

These relations contain the information to generate the tripartite networks to be analysed. They will be added in a folder "networks/*\_networks" where * is the name of the domains classification type (for example, superfamilies or FunFams), and the filename will be "network\_", followed by the annotation name. These networks contain two different types of layers: domain-protein and protein-FunSys, connecting domains with FunSys through protein nodes. These networks are generated with merge_pairs.rb.
*Note: add_protein_functional_families.rb output can be modified in the code if only a set of annotations and protein domains classification is provided.*

``` 
merge_pairs.rb -i domProtFunSysRelations -k domainTypeID -o network -n numberOfFiles -m minNumConns
```
Where: 
* -i domProtFunSysRelations = Input file with domains-proteins-FunSys relationships.
* -k domainTypeID = Domain type identifier. For example, if FunFams from CATH used, 'ff' identifier must be provided.
* -o tripartiteNetwork = Tripartite network output file.
* -n numberOfFiles = Number of files to output (if k-cross performed).
* -m minNumConns = Minimum number of proteins supporting a relation between domain and FunSys.

This *tripartiteNetwork* will contain the relations to be analysed by NetAnalyzer. This tool will return the list of domains associated with FunSys and the association values calculated depending on the index selected.

```
NetAnalyzer.rb -i tripartiteNetwork -l layers -m assocMeth -a assocVals
```
Where:
* -i tripartiteNetwork = Tripartite network to calculate associations between their nodes.
* -l layers = Layer names of the network and their identifier. For example, a tripartite network with domains classified in FunFams (ff), annotations in GO terms and proteins will be set as: "domains,ff;annotations,GO:;protID,[A-Za-z0-9]". Please separate between names and identifiers with commas, and between layers with semicolons.
* -m assocMeth = Association method (Jaccard, Pearson correlation coefficient or Hypergeometric index, amongst others). Set only one method by execution. For more information about association methods available in NetAnalyzer, please go to its documentation. 
* -U projection = Layers projection to calculate associations. In this case, we use proteins as nexus to calculate the association between domains and FunSys in common. For this, it should be established as "annotations,domains;protID". Please separate the nexus node (proteins) with semicolon, and nodes to associate (domains with FunSys) with comma.
* -a assocVals = Associations output file.

### STEP 3. Training DomFun with domain-FunSys associations and predicting proteins function ###
DomFun must be trained with associations between protein domains and FunSys to perform predictions. Depending on the association index set up in NetAnalyzer, these values are standardized or not before training DomFun. For example, in the case of using the Hypergeometric index, it is not necessary to perform this standardization. However, for the rest of methods is mandatory. For this, use standardize_scores.R.

*Note: please OMIT this step if hypergeometric index was used to calculate association values (they are directly transformed into P-values.*

```
standardize_scores.R -d assocVals -e threshold -o assocValsStd -s col > outZscore
```

Where
* -d assocVals = Domain-FunSys association values file to standardize.
* -e threshold = Association value to use as a threshold.
* -o assocValsStd = Association values standardized output file.
* -s col > outZscore = Z-score to use as a filter (if necessary).

DomFun will be trained with associations (standardized or not) to perform functional annotations predictions from a list of proteins (UniProt identifiers).

```
domains_to_function_predictor.rb -a assocVals -u -f cathDomainsFile -p proteinIDsList -c domainTypeID -T assocValThresh -i combMeth -t combScorThresh -o predictionResults
```
Where:
* -a assocVals = File with associations to train DomFun.
* -u multipleProts = Set if multiple proteins are given to predict. 
* -f cathDomainsFile = Domains file.
* -p proteinIDsList = List of protein identifiers.
* -c domainTypeID = Domain type identifier. For example, if FunFams from CATH used, 'ff' identifier must be provided.
* -T assocValThresh = Threshold to filter association values.
* -i combMeth = Method to combine association values. Please select between "fisher" and "stouffer".
* -t combScorThresh = Threshold to filter scores once combined.
* -o predictionResults = Output file with prediction results.

DomFun predictions output file (predictionResults) will contain a table with proteins predicted, their domains found, with which FunSys these domains where associated and the combined score calculated for the prediction.

The rest of scripts included in `bin/` folder were used to perform data transformation, validate results and plot graphs. They are described as follows:
* `association_metrics_average.rb`: Calculate average between different association files for results validation.
* `generate_CAFA2_dataset.rb`: Parse data from CAFA 2 validation system.
* `generate_CAFA2_tripartite_network.rb`: Generate tripartite networks from CAFA 2 dataset.
* `generate_cafa_control.rb`: Generate proteins control file from CAFA 2.
* `get_kegg_pathways.R`: Download KEGG pathways identifiers.
* `lines.R`: Program to plot different data distributions.
* `normalize_combined_scores.rb`: Script to normalize combined scores.
* `prepare_cafa_network.rb`: Generate a tripartite network with GO associations for CAFA 2 testing.
* `translate_kegg_genes2pathways.rb`: Use to translate KEGG identifiers from geneIDs to KEGG pathway IDs (from UniProt downloads).
* `validate_ProtFunSys_predictions.rb`: System to validate predictions when training DomFun with UniProt proteins.

## Development

After checking out the repo, run `bin/setup` to install dependencies. Then, run `rake spec` to run the tests. You can also run `bin/console` for an interactive prompt that will allow you to experiment.

To install this gem onto your local machine, run `bundle exec rake install`. To release a new version, update the version number in `version.rb`, and then run `bundle exec rake release`, which will create a git tag for the version, push git commits and tags, and push the `.gem` file to [rubygems.org](https://rubygems.org).

## Contributing

Bug reports and pull requests are welcome on BitBucket at https://bitbucket.org/elenarojano/domfun

## License

The gem is available as open source under the terms of the [MIT License](https://opensource.org/licenses/MIT).

