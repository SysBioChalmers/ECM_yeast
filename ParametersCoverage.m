%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Kcats,KMs] = ParametersCoverage(model_data,org_name)
%
% Kinetic Parameters coverage estimation (Kcat and Km values for substrates
% and products)for all of the enzymatic reactions  with an associated EC 
% number in model_data.model. 
%
% Ivan Domenzain.   Last edited: 2018-02-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ks,Kp] = ParametersCoverage(model_data,org_name,parameter)

    current       = pwd;
    EC_numbers    = model_data.EC_numbers;
    [nR, x]       = size(EC_numbers);
    Totalcounter  = 0;
    Kcatcounter   = 0;
    KScounter     = 0;
    KPcounter     = 0;
    model         = model_data.model;
   %Extract Kcat values, SA*Mw and KMs from the BRENDA files
    [data_cell] = extractBRENDAdata(parameter);
   %Creates a Structure with KEGG codes for organisms, names and phylogenetic 
   %distance matrix and extract the organism index in the KEGG struct
    cd ([current '/Databases'])
    load('TaxonomicDistance_KEGG.mat')
    org_index      = find_inKEGG(org_name,phylDistStruct.names);
    cd (current)
    KmTriplets = 0; KmNarrow = 0;
    KpTriplets = 0; KpNarrow = 0;
    
   %Create .txt files for saving the results         
    cd dataCoverage
    fileName_subs  = [parameter '_subs.txt'];
    fileID_subs    = fopen(fileName_subs,'w');
    fileName_prods = [parameter '_prods.txt'];
    fileID_prods   = fopen(fileName_prods,'w');
       
    %for each of the model reactions
    Kp = []; Ks = [];
    for i=1:nR
        subsIndx       = find(model.S(:,i)<0);
        prodsIndx      = find(model.S(:,i)>0);
        substrates     = model.metNames(subsIndx);
        nS             = length(substrates);
        products       = model.metNames(prodsIndx);
        nP             = length(products);
        nonEmpty       = find(~cellfun(@isempty, EC_numbers(i,:)));
        rxnECs         = unique(EC_numbers(i,nonEmpty));
        nE             = length(rxnECs);
       %For each of the possible isoenzymes for the i-th reaction
        for j=1:nE
            ECs     = strsplit(rxnECs{j},' ');         
           %For each of the enzymatic subunits (in case of complexes)
            Km_list = []; Kcat = [];
           %look for the parameters reported for each substrate in the rxn
            for k=1:nS
                cd ..
                [val,origin] = getKvalforMets(substrates(k),ECs,data_cell,...
                               phylDistStruct,org_index,org_name,parameter);
                
                str = ['rxn: ' num2str(i) ',[' rxnECs{j} '],' substrates{k} ',' num2str(val) ',' origin '\n'];
                if ~isnan(val)
                    Ks = [Ks; val];
                end
                cd dataCoverage
                fprintf(fileID_subs,str);
                disp(str)
            end
           %look for the parameters reported for each product in the rxn
            for k=1:nP
                cd ..
                [val,origin] = getKvalforMets(products(k),ECs,data_cell,...
                               phylDistStruct,org_index,org_name,parameter);
                
                str = ['rxn: ' num2str(i) ',[' rxnECs{j} '],' products{k} ',' num2str(val) ',' origin '\n'];
                if ~isnan(val)
                    Kp = [Kp; val];
                end
                cd dataCoverage
                fprintf(fileID_prods,str);
                disp(str)
            end
             
        end                
        %Paremeters needed by the model
        Kcatcounter   = Kcatcounter + nE;
        KScounter     = KScounter + nE*nS;
        KPcounter     = KPcounter + nE*nP;
        Totalcounter  = Totalcounter + nE + nE*(nS+nP);
        
    end
    cd (current)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val,origin] = getKvalforMets(met,ECs,data_cell,DistStruct,indx,...
                                                        org_name,parameter)
    options = {org_name;'closest'};
    for l=1:length(ECs)
            val = []; ii = 1;
            while isempty(val) && ii<=2
                [val,KmTriplets,KmNarrow] = extractKval(ECs(l),data_cell,...
                                         met,indx, DistStruct,options(ii));
                ii=ii+1;
            end

            if ~isempty(val)
                %Take the minimum Kcat value to avoid underpredicted 
                %enzyme costs
                if strcmpi(parameter,'kcat')
                    val = max(val);
                elseif strcmpi(parameter,'km')
                %Take the max Km values to avoid underpredicted 
                %enzyme costs
                    val = min(val);
                end
            else
                val = NaN;
            end
    end
    origin = options{ii-1};
end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Load BRENDA data:
 function [data_cell] = extractBRENDAdata(parameter)
    cd Databases
    if strcmpi(parameter,'kcat')
        %Get all the Kcat values (maximum values for the organism/substrate/EC#
        %triplets distribution) for non-mutant enzymes.
        disp('Extracting Kcat values')
        fID           = fopen('max_KCAT.txt');
        data_cell     = textscan(fID,'%s %s %s %f  %s','delimiter','\t');
        data_cell{4}  = data_cell{4};  
        fclose(fID);
    elseif strcmpi(parameter,'km')
        %Get all the KM values (minimum values for the organism/substrate/EC#
        %triplets distribution) for non-mutant enzymes.
        disp('Extractring KM values')
        fID         = fopen('min_KM.txt');
        data_cell     = textscan(fID,'%s %s %s %f  %s','delimiter','\t');
        data_cell{4}  = data_cell{4}*1000;   %[microM] -> [miliM]
        fclose(fID);
    end
    %Split string for each organism in the BRENDA data 
    %{name, taxonomy, KEGG code}
    data_cell{3}   = cellfun(@stringSplit, data_cell{3});
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
        %If no KEGG code was found then it will look for first match on
        %comparing genus names
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Kval,Triplets,NarrowDist] = extractKval(EC,Kcell,met,OrgIndx,...
                                                  DistStruct,org_name)
    Kval = []; NarrowDist = 0;Triplets = 0;                                           
   % Extract indexes for the EC# appearences for any met or
   % organism)
    EC_indexes = extract_indexes(EC,Kcell,'*',false,'any', OrgIndx,...
                                 DistStruct,false);
    if ~isempty(EC_indexes)
    	values   = Kcell{4}(EC_indexes);
        wideness = abs(log(median(values))/...
                       log(max(values)/min(values)));
       %If the Kvalues present a narrow distribution then
       %the mean value is taken
        if wideness >= 1
        	%Kval       = [Kval; mean(values)];
            Kval       = mean(values);
            NarrowDist = NarrowDist + 1;
        else
       %Search for the triplet [EC#/Substrate/Organism] on the 
       %KM_cell data
            EC_indexes = extract_indexes(EC,Kcell,met,true,org_name, ...
                                                 OrgIndx,DistStruct,false);

            if ~isempty(EC_indexes)
                %Kval     = [Kval; min(Kcell{4}(EC_indexes))];
                Kval     = Kcell{4}(EC_indexes);
                Triplets = Triplets + 1;
            end
        end

    else
        %disp('EC number not found in data')            
    end        
end 
