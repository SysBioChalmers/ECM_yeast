function model  = createKEGGmodel(model,yeastGEM)
%createKEGGmodel
%  
% Converts a MATLAB model structure into a .txt KEGG model in rxns are
% explicitly expressed with their formulas and matched to KEGG rxn IDs
% retrieved from yeast8, if available.
%
%     model   Initial MATLAB model structure (subset of yeastGEM)
%     yeast8  MATLAB structure of the latest yeastGEM
%
% Ivan Domenzain.   Last edited: 2019-10-11
%

model = getKEGGmets(model,yeastGEM);
model = mapRxnInfo(model,yeastGEM);
model = buildKEGGrxns(model);
%Get the rxn indexes with gene(s) association
Indxs           = ~cellfun(@isempty,model.grRules);
ReactionFormula = model.KEGGrxns(Indxs);
ID              = model.rxns(Indxs);
Name            = model.rxnNames(Indxs);
KEGGIDs         = model.rxnKEGGID(Indxs);
%!ID	!Name	!ReactionFormula	!Identifiers:kegg.reaction		
KEGGtable = table(ID,Name,ReactionFormula,KEGGIDs);
writetable(KEGGtable,'../../models/KEGG_model_Scaffold.txt','Delimiter','\t')
end