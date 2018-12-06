function model = addNewFields(model,yeast)
%Creates or modifies the fields grRules, rxnKEGGID, rxnGeneMat consitently
%with the yeast8 model
model = mapRxnInfo(model,yeast);
%Creates the field metKEGGID (KEGG ids for every metabolite in the network
model = getKEGGmets(model,yeast);
%Adds the field KEGGrxns which includes the rxn equations with KEGG IDs for
%the compounds
model = buildKEGGrxns(model);
end
