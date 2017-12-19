%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Kcats,] = ParametersCoverage(model_data,modelECs)
%
% Kinetic Parameters coverage estimation (Kcats and Km values for substrates
% and products)for all of the enzymatic reactions  with an associated EC 
% number in model_data.model. 
%
% Ivan Domenzain.   Last edited: 2017-12-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [Kcats,] = ParametersCoverage(model_data,modelECs)

    EC_numbers    = model_data.EC_numbers;
    model         = model_data.model;
    [nR, x]       = size(EC_numbers);
    Totalcounter  = 0;
    Kcatcounter   = 0;
    KScounter     = 0;
    KPcounter     = 0;
   %Creates a Structure with KEGG codes for organisms, names and phylogenetic 
   %distance matrix and extract the organism index in the KEGG struct
    phylDistStruct =  KEGG_struct(true);
    org_index      = find_inKEGG(org_name,phylDistStruct.names);
    %Extract Kcat values, SA*Mw and KMs from the BRENDA files
    [Kcat_cell, SA_cel, KM_cell] = extractBRENDAdata;
    for i=1:nR
        subsIndx      = find(model.S(:,i)<0);
        prodsIndx     = find(model.S(:,i)>0);
        substrates    = model.metNames(subsIndx);
        products      = model.metNames(prodsIndx);
        nonEmpty      = find(~cellfun(@isempty, EC_numbers(i,:)));
        rxnECs        = unique(EC_numbers(i,nonEmpty));
        Kcatcounter   = Kcatcounter + length(rxnECs);
        KScounter     = KScounter + length(rxnECs)*(length(substrates));
        KPcounter     = KPcounter + length(rxnECs)*(length(products));
        Totalcounter  = Totalcounter + length(rxnECs) + ...
                        length(rxnECs)*(length(substrates) + length(products));


        % Next step! to read data from 

    end
 %end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Load BRENDA data:
 function [Kcat_cell,SA_cell,KM_cell] = extractBRENDAdata
    cd Databases
    %Get all the Kcat values (maximum values for the organism/substrate/EC#
    %triplets distribution) for non-mutant enzymes.
    fID           = fopen('max_KCAT.txt');
    Kcat_cell     = textscan(fID,'%s %s %s %f  %s','delimiter','\t');
    Kcat_cell{4}  = Kcat_cell{4};   %[1/s] -> [1/h]
    fclose(fID);
    %Get all the KM values (minimum values for the organism/substrate/EC#
    %triplets distribution) for non-mutant enzymes.
    fID         = fopen('min_KM.txt');
    KM_cell     = textscan(fID,'%s %s %s %f  %s','delimiter','\t');
    KM_cell{4}  = KM_cell{4}*1000;   %[1/s] -> [1/h]
    fclose(fID);
    %Split string for each organism in the BRENDA data 
    %{name, taxonomy, KEGG code}
    Kcat_cell{3} = cellfun(@stringSplit, Kcat_cell{3});
    KM_cell{3}   = cellfun(@stringSplit, KM_cell{3});
    %Creates a cell with all the SA*Mw values (corresponding to the same EC#
    %and organism).
    %SA_cell      = SA_BRENDA('max_SA.txt','max_MW.txt');
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function string_cells = stringSplit(cell_array)
        string_cells = {strsplit(cell_array,'//')};
        string_cells = string_cells{1}(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function org_index = find_inKEGG(org_name,names)
    org_index      = find(strcmpi(org_name,names));
    if isempty(org_index)
        i=1;
        while isempty(org_index) && i<length(names)
            str = names{i};
            if strcmpi(org_name(1:strfind(org_name,' ')-1),...
                str(1:strfind(str,' ')-1))
                org_index = i;
            end
            i = i+1;
        end
        if isempty(org_index);org_index = '*';end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
