function model  = createKEGGmodel(model,yeast8)
cd Model_curation_scripts
model = getKEGGmets(model,yeast8);
model = mapRxnInfo(model,yeast8);
model = buildKEGGrxns(model);
%Get the rxn indexes with gene(s) association
Indxs  = ~cellfun(@isempty,model.grRules);
ReactionFormula   = model.KEGGrxns(Indxs);
ID      = model.rxns(Indxs);
Name    = model.rxnNames(Indxs);
KEGGIDs = model.rxnKEGGID(Indxs);
%!ID	!Name	!ReactionFormula	!Identifiers:kegg.reaction		
KEGGtable = table(ID,Name,ReactionFormula,KEGGIDs);
cd ../Models
writetable(KEGGtable,'KEGGmodelScaffold.txt','Delimiter','\t')
end