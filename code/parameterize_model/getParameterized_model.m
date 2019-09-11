function ECM_model =  getParameterized_model
%
% Ivan Domenzain.   Last edited: 2019-10-11

current = pwd;
%Specify the GECKO toolbox directory
repoPath_terminal = '/Users/ivand/Documents/GitHub/ECM_yeast';
yeastGEM          = cloneAndLoad('yeast-GEM','GECKO');
org_name          = 'saccharomyces cerevisiae';
org_code          = 'sce';
%% Preprocess reduced GEM  
model = load('../../models/model.mat');
model = model.model;
model = createKEGGmodel(model,yeastGEM);
model = addNewFields(model,yeastGEM);
save('../../models/reducedYeast.mat','model')
%Save SBML file
exportModel(model,'../../models/reducedYeast.xml',true)
%% Get parameters from BRENDA 
%Get EC numbers for each reaction based on the grRules field
cd ('../GECKO/geckomat/get_enzyme_data')
model_data = getEnzymeCodes(model);
cd (current)
%Get Kcats (forward and backward) for each reaction based on EC#
[KC_s, KC_p] = getKineticData(model_data,org_name,'catalytic rate constant');
%Get KM's (forward and backward) for each reaction based on EC#
[KM_s, KM_p] = getKineticData(model_data,org_name,'Michaelis constant');
cd ../visualization
plotCDF(KC_s,KC_p,'Kcat distributions','Kcat [1/s]','../../results/plots/collected_Kcats.fig')
plotCDF(KM_s,KM_p,'KM distributions','KM [mM]','../../results/plots/collected_KMs.fig')
%% Get Keq from Equilibrator  
%Use equilibrator (python) for determining Keq's for all of the enzymatic
%reactions
commandStr = ['python ' repoPath_terminal '/code/parameterize_model/getKeq.py'];
status     = system(commandStr);
if status~=0
     disp('getKeq.py failed to execute from MATLAB, try to run it separately')
     disp('in a terminal, once this is done, press any key to continue')
     pause
end
%% Get ECM model structure 
cd (current)
ECM_model = getECM_model(model);
cd ..
rmdir('GECKO', 's')
rmdir('yeast-GEM', 's')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = cloneAndLoad(repo1,repo2)
cd ..
git(['clone https://github.com/SysBioChalmers/' repo1])
git(['clone https://github.com/SysBioChalmers/' repo2 '.git'])
load([repo1 '/ModelFiles/mat/yeastGEM.mat'])
%model = ravenCobraWrapper(model);
cd parameterize_model
clc
end