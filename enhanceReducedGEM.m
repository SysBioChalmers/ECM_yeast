%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model_data,kcats,KMs] = enhanceReducedGEM(model,org_name,KEGGcode)
%
% Ivan Domenzain.   Last edited: 2018-10-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [model_data,kcats,KMs] = enhanceGEM(model,org_name,org_code)
 current       = pwd; 
 %Specify the GECKO toolbox directory
 repoPath_terminal = '/Users/ivand/Documents/GitHub/ECM_yeast';
 GECKO_path    = '/Users/ivand/Documents/GitHub/GECKO';
 yeast8_path   = '/Users/ivand/Documents/GitHub/yeast-GEM/ModelFiles/xml';
 org_name      = 'saccharomyces cerevisiae';
 org_code      = 'sce';
 %%%%%%%%%%%%%%%%%%%%%%% Preprocess reduced GEM  %%%%%%%%%%%%%%%%%%%%%%%%%%
 %Commented because the missing genes were also manually added to the model
 cd (yeast8_path)
 yeast8 = readCbModel('yeastGEM.xml');
 cd ([current '/Models'])
 model = importExcelModel('model.xlsx');
 cd ..
 model = createKEGGmodel(model,yeast8);
 %model = correctMetNames(model);
 cd ../Model_curation_scripts
 model = addNewFields(model,yeast8);
 % Get EC numbers for each reaction based on the grRules field
 cd ([GECKO_path '/geckomat/get_enzyme_data'])
 model_data = getEnzymeCodes(model);
 cd ([current '/Models'])
 save('reducedModel.mat','model','model_data') 
 %Save SBML file
 exportModel(model,'reducedYeast.xml',true)
 %%%%%%%%%%%%%%%%%%%% Get parameters from BRENDA %%%%%%%%%%%%%%%%%%%%%%%%%%
 cd (current)
 %Get all the unique EC# matched to the model and then analyse the distribution
 %of KM and Kcat values reported for each of these. EC numbers with more
 %than one value and narrow distributions (log(median)/log(max/min)>= 1 are
 %showed in the cell 'NarrowDists'
 %[NarrowKcats,KcatStats] = BRENDA_analysis(model_data,'Kcat');
 %[NarrowKMs,KMStats]     = BRENDA_analysis(model_data,'KM');
 %Reversibility consistency check between the LB and rev fields in the model
 %inconsistencies          = reversibilityConsistencyCheck(model);
 [Ks, Kp] = getKineticData(model_data,org_name,'catalytic rate constant');
 %plotCDF(Ks,Kp,'Kcat distributions','Kcat [1/s]')
 [Ks, Kp] = getKineticData(model_data,org_name,'Michaelis constant');
 %plotCDF(Ks,Kp,'KM substrates','KM [mM]')

%%%%%%%%%%%%%%%%%%%% Get Keq from Equilibrator  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use equilibrator (python) for determining Keq's for all of the enzymatic
%reactions
commandStr = ['python ' repoPath_terminal '/getKeq.py'];
status = system(commandStr);
if status~=0
    disp('getKeq.py failed to execute from MATLAB, try to run it separately')
    disp('in a terminal, once this is done, press any key to continue')
    pause
end
%%%%%%%%%%%%%%%%%%%%%%%% Get ECM model structure %%%%%%%%%%%%%%%%%%%%%%%%%%
ECM_model = getECM_model(model);
cd (current)
%end

