%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecModel,model_data,kcats] = enhanceGEM(model,toolbox,name)
%
% Benjam?n J. S?nchez. Last edited: 2017-04-12
% Ivan Domenzain.   Last edited: 2017-11-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [ecModel,model_data,kcats] = enhanceGEM(model,toolbox,org_name,org_code)
 GECKO_path    = '/Users/ivand/Documents/GECKO-IVAN';
 Protdatabase  = 'sce_ProtDatabase.mat';
 org_name      = 'saccharomyces cerevisiae';
 org_code      = 'sce';
 gR_exp        = 0.41;   % Batch growth on minimal glucose media [g/gDw hr]
  toolbox      = 'COBRA';
%  format short e
%  if strcmp(toolbox,'COBRA')
%     initCobraToolbox;
%  end
  cd ([GECKO_path '/Models'])
% % %sce_model = importModel('yeast_7.00_cobra.xml',true,true)
%  %Update yeast 7 model with all recommended changes:
  cd ../Matlab_Module/get_enzyme_data
   %model = modelCorrections(model);
% % %Add some RAVEN fields for easier visualization later on:
%  if strcmpi(org_name,'saccharomyces cerevisiae')
%      model = standardizeModel(model,toolbox);
%  end 
% %Retrieve kcats & MWs for each rxn in model:
model_data               = getEnzymeCodes(model,Protdatabase,GECKO_path);
[UniqueECs, NarrowDists] = BRENDA_analysis(model_data);

%KMs        = matchKcats_AVLANT(model_data, org_name);

%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
