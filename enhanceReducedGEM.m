%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model_data,kcats,KMs] = enhanceReducedGEM(model,org_name,KEGGcode,...
%                                                             Protdatabase)
%
% Ivan Domenzain.   Last edited: 2017-12-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [model_data,kcats,KMs] = enhanceGEM(model,toolbox,org_name,org_code)
 current       = pwd; 
 %Specify the GECKO toolbox directory
 GECKO_path    = '/Users/ivand/Documents/GECKO-IVAN';
 Protdatabase  = 'sce_ProtDatabase.mat';
 %org_name      = 'saccharomyces cerevisiae';
 %org_code      = 'sce';
 toolbox       = 'COBRA'; 
 cd ([GECKO_path '/Matlab_Module/get_enzyme_data'])
 % Get EC numbers for each reaction based on the grRules field
 model_data               = getEnzymeCodes(model,Protdatabase,GECKO_path);
 cd (current)
 %Get the unique EC# matched to the model and then analyse the distribution
 %of KM and Kcat values reported for each of these. EC numbers with more
 %than one value and narrow distributions (log(median)/log(max/min)>= 1 are
 %showed in the cell 'NarrowDists'
 [NarrowKcats,KcatStats] = BRENDA_analysis(model_data,'Kcat');
 [NarrowKMs,KMStats]     = BRENDA_analysis(model_data,'KM');
 %Reversibility consistency check between the LB and rev fields in the model
 inconsistencies          = reversibilityConsistencyCheck(model);
 
 ParametersCoverage(model_data,org_name,'kcat')
 ParametersCoverage(model_data,org_name,'km')

%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
