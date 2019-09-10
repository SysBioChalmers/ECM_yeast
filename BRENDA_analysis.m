%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ECnumbers, NarrowDists] = BRENDA_Analysis(model_data)
%
% Gets a txt file generated with the script findMaxKvalues 
% (GECKO/Python_Module) and analyses the distribution of Kvalues 
% (Kcats,KM, SA, Mw) for each of the EC numbers in the model_data.
%
% Ivan Domenzain.   Last edited: 2017-12-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Narrow, stats] = BRENDA_analysis(model_data, parameter)
    % Extract all the unique EC numbers included in the model
    ECnumbers = extractECnumbers(model_data.EC_numbers);
    cd Databases
    
    if strcmpi(parameter,'KCAT')
        file             = 'max_KCAT.txt';
        conversionFactor = 1;    %[1/s] ->[1/s]
    elseif strcmpi(parameter,'KM')
        file             = 'min_KM.txt';
        conversionFactor = 1000; %[microM] ->[miliM]
    end
    
    fID     = fopen(file);
    data    = textscan(fID,'%s %s %s %f  %s','delimiter','\t');
    data{4} = data{4}*conversionFactor;
    fclose(fID);
    %Split string for each organism in the BRENDA data {name,taxonomy,KEGG code}
    data{3} = cellfun(@stringSplit, data{3});
    %The parameters matching criteria can be flexibilized for the EC numbers
    %with narrow distributions
    str            = [parameter ' distributions wideness'];
    [Narrow,stats] = analyseDistributions(ECnumbers,data{1},data{4},str);
    cd ..
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wideDists, stats] = analyseDistributions(ECs,EC_cell,Kvalues,str)
    wideness  = [];
    wideDists = []; 
    for i=1:length(ECs)
         indx      = find(strcmpi(ECs(i),EC_cell));
         if ~isempty(indx)
             EC_KMdist = Kvalues(indx);
             if length(indx) == 1
                 metric = 0;
             else
                 metric = log10(median(EC_KMdist))/log10(max(EC_KMdist)/...
                                                           min(EC_KMdist));
                 wideness  = [wideness; metric];
             end
             if metric >= 1
                 wideDists = [wideDists; ECs(i)]; 
             end
         end

    end
    figure
    [y, stats]=cdfplot(wideness);
    title(str)
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


