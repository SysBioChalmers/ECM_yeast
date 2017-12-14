%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BrendaDataAnalysis(model_data)
%
% Gets a txt file generated with the findMaxKvalues (GECKO/Python_Module)
% and analyses the distribution of Kvalues (Kcats,KM, SA, Mw) for each of
% the EC numbers in the model_data.
%
% Ivan Domenzain.   Last edited: 2017-12-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract all the unique EC numbers included in the model
ECnumbers = extractECnumbers(model_data.EC_numbers);
cd Databases
%kcat_file = 'max_KCAT.txt';
km_file   = 'min_KM.txt';
%SA_file    = 'max_SA.txt';
%MW_file    = 'max_MW.txt';
fID       = fopen(km_file);
KM        = textscan(fID,'%s %s %s %f  %s','delimiter','\t');
KM{4}     = KM{4}*1000;   %[mM] -> [microM]
fclose(fID);
%Split string for each organism in the BRENDA data {name,taxonomy,KEGG code}
KM{3}     = cellfun(@stringSplit, KM{3});
% The parameters matching criteria can be flexibilized for the EC numbers
% with narrow distributions
[NarrowDists, stats] = getWidenessDistributions(ECnumbers,KM{1},KM{4});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wideDists, stats] = getWidenessDistributions(ECs,EC_cell,Kvalues)
    wideness  = [];
    wideDists = []; 
    for i=1:length(ECs)
         indx      = find(strcmpi(ECs(i),EC_cell));
         if ~isempty(indx)
             EC_KMdist = Kvalues(indx);
             metric    = log10(mean(EC_KMdist))/log10(max(EC_KMdist)/...
                                                          min(EC_KMdist));
             wideness  = [wideness; metric];
             if metric >= 1
                 wideDists = [wideDists; ECs(i)]; 
             end
         end

    end
    [y, stats]=cdfplot(wideness);
    ylabel('Cumulative distribution','FontSize',30,'FontWeight','bold');
    xlabel('log10(mean)/log10(max/min)','FontSize',30,'FontWeight','bold');
    hold on

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function string_cells = stringSplit(cell_array)
        string_cells = {strsplit(cell_array,'//')};
        string_cells = string_cells{1}(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ECnumbers = extractECnumbers(EC_numbers)

    [M,N]     = size(EC_numbers);
    ECnumbers =  [];
    for i=1:N
        for j=1:M
            if ~isempty(EC_numbers{j,i})
                datum = strsplit(EC_numbers{j,i},' ');
                for k=1:length(datum)
                    if ~ismember(datum{k},ECnumbers)
                      ECnumbers = [ECnumbers; datum(k)];  
                    end
                end
            end
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


