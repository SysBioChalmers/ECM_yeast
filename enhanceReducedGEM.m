%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model_data,kcats,KMs] = enhanceReducedGEM(model,org_name,KEGGcode)
%
% Ivan Domenzain.   Last edited: 2018-04-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [model_data,kcats,KMs] = enhanceGEM(model,org_name,org_code)
 current       = pwd; 
 %Specify the GECKO toolbox directory
 repoPath_terminal = '/Users/ivand/Documents/GitHub/ECM_yeast';
 GECKO_path    = '/Users/ivand/Documents/GitHub/GECKO';
 yeast8_path   = '/Users/ivand/Documents/GitHub/yeast-GEM/ModelFiles/xml';
 org_name      = 'saccharomyces cerevisiae';
 org_code      = 'sce';
 %%%%%%%%%%%%%%%%%%%%%%%%%% Model modifications  %%%%%%%%%%%%%%%%%%%%%%%%%%
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
 %%%%%%%%%%%%%%%%%%%% Parameter distributions analysis  %%%%%%%%%%%%%%%%%%%
 cd (current)
 %Get all the unique EC# matched to the model and then analyse the distribution
 %of KM and Kcat values reported for each of these. EC numbers with more
 %than one value and narrow distributions (log(median)/log(max/min)>= 1 are
 %showed in the cell 'NarrowDists'
 %[NarrowKcats,KcatStats] = BRENDA_analysis(model_data,'Kcat');
 %[NarrowKMs,KMStats]     = BRENDA_analysis(model_data,'KM');
 %Reversibility consistency check between the LB and rev fields in the model
 %inconsistencies          = reversibilityConsistencyCheck(model);
 [Ks, Kp] = ParametersCoverage(model_data,org_name,'catalytic rate constant');
 %plotCDF(Ks,Kp,'Kcat distributions','Kcat [1/s]')

 [Ks, Kp] = ParametersCoverage(model_data,org_name,'Michaelis constant');
 %plotCDF(Ks,Kp,'KM substrates','KM [mM]')

%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use equilibrator (python) for determining Keq's for all of the enzymatic
%reactions
commandStr = ['python ' repoPath_terminal '/build_ECM_modelFile.py'];
 system(commandStr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotCDF(distribution1,distribution2,name,xAxis)
    figure
    [y, stats]=cdfplot(distribution1);
    hold on
    [y, stats]=cdfplot(distribution2);
    hold on
    title(name)
    ylabel('Cumulative distribution','FontSize',30,'FontWeight','bold');
    xlabel(xAxis,'FontSize',30,'FontWeight','bold');
    legend('Substrates','Products')
end