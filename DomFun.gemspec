
lib = File.expand_path("../lib", __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require "DomFun/version"

Gem::Specification.new do |spec|
  spec.name          = "DomFun"
  spec.version       = DomFun::VERSION
  spec.authors       = ["Elena Rojano, Pedro Seoane"]
  spec.email         = ["elenarojano@uma.es, seoanezonjic@hotmail.com"]

  spec.summary       = %q{Tool to predict protein functions using domains-FunSys associations.}
  spec.description   = %q{Proteins function predictor trained with associations between protein domains (classified in CATH) and functional systems (GO, KEGG, Reactome).}
  spec.homepage      = "https://bitbucket.org/elenarojano/domfun"
  spec.license       = "MIT"

  # Prevent pushing this gem to RubyGems.org. To allow pushes either set the 'allowed_push_host'
  # to allow pushing to a single host or delete this section to allow pushing to any host.
  # if spec.respond_to?(:metadata)
  #   spec.metadata["allowed_push_host"] = "TODO: Set to 'http://mygemserver.com'"

  #   spec.metadata["homepage_uri"] = spec.homepage
  #   spec.metadata["source_code_uri"] = "TODO: Put your gem's public repo URL here."
  #   spec.metadata["changelog_uri"] = "TODO: Put your gem's CHANGELOG.md URL here."
  # else
  #   raise "RubyGems 2.0 or newer is required to protect against " \
  #     "public gem pushes."
  # end

  # Specify which files should be added to the gem when it is released.
  # The `git ls-files -z` loads the files in the RubyGem that have been added into git.
  spec.files         = Dir.chdir(File.expand_path('..', __FILE__)) do
    `git ls-files -z`.split("\x0").reject { |f| f.match(%r{^(test|spec|features)/}) }
  end
  spec.bindir        = "exe"
  spec.executables   = spec.files.grep(%r{^exe/}) { |f| File.basename(f) }
  spec.require_paths = ["lib"]

  spec.add_development_dependency "bundler", "~> 2.0"
  spec.add_development_dependency "rake", "~> 12.3.1"
  spec.add_development_dependency "rspec", "~> 3.0"
  
  spec.add_dependency "NetAnalyzer", "~> 0.1.5"
  spec.add_dependency "statistics2", ">= 0.54"
  spec.add_dependency "terminal-table", ">= 2.0.0"
  spec.add_dependency "report_html", ">= 0.4.3"
  spec.add_dependency "parallel"

end
