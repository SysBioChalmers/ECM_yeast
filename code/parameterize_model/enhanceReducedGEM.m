function ECM_model =  enhanceReducedGEM
%enhanceReducedGEM
%
% Ivan Domenzain.   Last edited: 2019-10-11

current = pwd;
%Specify the GECKO toolbox directory
repoPath_terminal = '/Users/ivand/Documents/GitHub/ECM_yeast';
GECKO_path    = '/Users/ivand/Documents/GitHub/GECKO';
yeast8_path   = '/Users/ivand/Documents/GitHub/yeast-GEM/ModelFiles/mat';
org_name      = 'saccharomyces cerevisiae';
org_code      = 'sce';
%%%%%%%%%%%%%%%%%%%%%%% Preprocess reduced GEM  %%%%%%%%%%%%%%%%%%%%%%%%%%
model = load('../../models/model.mat');
model = model.model;
yeast8 = load([yeast8_path '/yeastGEM.mat']);
yeast8 = yeast8.model;
model = createKEGGmodel(model,yeast8);
model = addNewFields(model,yeast8);
save('../../models/reducedModel.mat','model')
%Save SBML file
exportModel(model,'../../models/reducedYeast.xml',true)
%%%%%%%%%%%%%%%%%%%% Get parameters from BRENDA %%%%%%%%%%%%%%%%%%%%%%%%%%
% Get EC numbers for each reaction based on the grRules field
cd ([GECKO_path '/geckomat/get_enzyme_data'])
model_data = getEnzymeCodes(model);
cd (current)
%Get Kcats (forward and backward) for each reaction based on EC#
[Ks, Kp] = getKineticData(model_data,org_name,'catalytic rate constant');
cd ../visualization
plotCDF(Ks,Kp,'Kcat distributions','Kcat [1/s]')
%Get KM's (forward and backward) for each reaction based on EC#
cd ../parameterize_model
[Ks, Kp] = getKineticData(model_data,org_name,'Michaelis constant');
cd ../visualization
plotCDF(Ks,Kp,'KM substrates','KM [mM]')
%%%%%%%%%%%%%%%%%%%% Get Keq from Equilibrator  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use equilibrator (python) for determining Keq's for all of the enzymatic
%reactions
cd ../parameterize_model
commandStr = 'python getKeq.py';
status = system(commandStr);
if status~=0
    disp('getKeq.py failed to execute from MATLAB, try to run it separately')
    disp('in a terminal, once this is done, press any key to continue')
    pause
end
%%%%%%%%%%%%%%%%%%%%%%%% Get ECM model structure %%%%%%%%%%%%%%%%%%%%%%%%%%
ECM_model = getECM_model(model);
cd (current)
end

