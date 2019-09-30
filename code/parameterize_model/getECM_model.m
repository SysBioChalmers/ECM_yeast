function ECM_model = getECM_model(model)
%
% getECM_model
%
% Get a GEM matlab structure and transforms it into an ECM metabolic
% network.
%
% usage: ECM_model = getECM_model(model)
%
% Ivan Domenzain.   Last edited: 2019-09-11
%

%ECM operates with exclusively with eznymatic reactions
enzRxns = find(~cellfun(@isempty,model.grRules));
nR      = length(enzRxns);
nM      = length(model.mets);
S       = full(model.S(:,enzRxns));
rev     = ones(nR,1);
mets    = model.mets;
rxns    = model.rxns(enzRxns);
%Indicate external metabolites
ext       = [];
kinetics  = true;
ECM_model = network_construct(S,rev,ext,mets,rxns,kinetics);
%Add missing fields
ECM_model.genes                   = model.genes(enzRxns);
ECM_model.reaction_KEGGID         = model.rxnKEGGID(enzRxns);
ECM_model.reaction_NameForPlots   = model.rxns(enzRxns);
ECM_model.metabolite_names        = model.metNames;
ECM_model.metabolite_NameForPlots = model.mets;
ECM_model.metabolite_KEGGID       = model.metKEGGID;
%Get and save an SBML file for the ECM_model (including kinetic
%expressions)
SBMLmodel = network_sbml_export(ECM_model,'','reducedYeast_ECM','../../models/reducedYeast_ECM.xml');
%Create network structure
options = struct;
options.filename                      = '../../models/reducedYeast_ECM.tsv';
options.modular_rate_law_table        = enzRxns;
options.modular_rate_law_kinetics     = enzRxns;
options.modular_rate_law_parameter_id = enzRxns;
options.save_in_one_file              = 1;
options.document_name                 = 'reducedYeast_ECM';
%Get and save an SBtab version of the ECM model
sbtab_document = network_to_sbtab(ECM_model, options);
%addParametersToNetwork(options.document_name)
end
