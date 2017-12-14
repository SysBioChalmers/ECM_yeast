%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecModel,model_data,kcats] = enhanceGEM(model,toolbox,name)
%
% Ivan Domenzain.   Last edited: 2017-12-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [ecModel,model_data,kcats] = enhanceGEM(model,toolbox,org_name,org_code)
 GECKO_path    = '/Users/ivand/Documents/GECKO-IVAN';
 Protdatabase  = 'sce_ProtDatabase.mat';
 org_name      = 'saccharomyces cerevisiae';
 org_code      = 'sce';
 cd ../Matlab_Module/get_enzyme_data
% Get EC numbers for each reaction based on the grRules field
model_data               = getEnzymeCodes(model,Protdatabase,GECKO_path);
[UniqueECs, NarrowDists] = BRENDA_analysis(model_data);
inconsistencies          = reversibilityConsistencyCheck(model);



%KMs        = matchKcats_AVLANT(model_data, org_name);

%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
