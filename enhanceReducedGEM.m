%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model_data,kcats,KMs] = enhanceReducedGEM(model,org_name,KEGGcode)
%
% Ivan Domenzain.   Last edited: 2018-02-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [model_data,kcats,KMs] = enhanceGEM(model,org_name,org_code)
 current       = pwd; 
 %Specify the GECKO toolbox directory
 GECKO_path    = '/Users/ivand/Documents/GitHub/GECKO';
 org_name      = 'saccharomyces cerevisiae';
 org_code      = 'sce';

 %model modifications
 model = AddMissingGenes(model); 
 model = correctMetNames(model);
 model = modelCorrections(model);
 % Get EC numbers for each reaction based on the grRules field
 cd ([GECKO_path '/Matlab_Module/get_enzyme_data'])
 model_data    = getEnzymeCodes(model);
 cd (current)
 cd Models
 save('reducedModel.mat','model','model_data') 
 cd ..
 %Get all the unique EC# matched to the model and then analyse the distribution
 %of KM and Kcat values reported for each of these. EC numbers with more
 %than one value and narrow distributions (log(median)/log(max/min)>= 1 are
 %showed in the cell 'NarrowDists'
 [NarrowKcats,KcatStats] = BRENDA_analysis(model_data,'Kcat');
 [NarrowKMs,KMStats]     = BRENDA_analysis(model_data,'KM');
 %Reversibility consistency check between the LB and rev fields in the model
 inconsistencies          = reversibilityConsistencyCheck(model);
 [Ks, Kp] = ParametersCoverage(model_data,org_name,'catalytic rate constant');
 plotCDF(Ks,Kp,'Kcat distributions','Kcat [1/s]')

 [Ks, Kp] = ParametersCoverage(model_data,org_name,'Michaelis constant');
 plotCDF(Ks,Kp,'KM substrates','KM [mM]')

%end
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