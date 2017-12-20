%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function EC_indexes = extract_indexes(EC,data_cell,metName,met_flag,...
%                                   organism, org_index,TaxDistStruct,wild)
%
% Extract the indexes of the entries in the BRENDA data that meet the 
% conditions specified by the search criteria
%
% Ivan Domenzain.   Last edited 2017-12-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EC_indexes = extract_indexes(EC,data_cell,metName,met_flag,...
                                    organism, org_index,TaxDistStruct,wild)
     
                     
    EC_cell   = data_cell{1}; 
    subs_cell = data_cell{2};
    orgs_cell = data_cell{3};
    
    EC_indexes = [];
    %In the case of wild cards present in the EC number
    if wild
      for j=1:length(EC_cell)
           if strfind(EC_cell{j},EC)==1
               EC_indexes = [EC_indexes,j];
           end
       end   
    else
       EC_indexes = transpose(find(strcmpi(EC,EC_cell)));  
    end
   %If met_flag=true then it will extract only the entry indexes for the 
   %specific metabolitein the EC subset from the BRENDA cell array
    if met_flag
        Subs_indexes = [];
        for l = 1:length(metName)
            if ~isempty(metName(l))
                Subs_indexes = horzcat(Subs_indexes,...
                    EC_indexes(strcmpi(metName(l),subs_cell(EC_indexes))));          
            end
        end
        EC_indexes = Subs_indexes;    
    end
    EC_orgs    = orgs_cell(EC_indexes);
 
   %If specific organism values are requested look for all the organism
   %repetitions on the subset BRENDA cell array(EC_indexes)
    if ~ismember(string(organism),{'closest';'any'}) 
        EC_indexes = EC_indexes(strcmpi(string(organism),EC_orgs));

   %If KEGG code was assigned to the organism (model) then it will look for   
   %the Kcat value for the closest organism
    elseif strcmpi(organism,'closest') && org_index~='*'
        KEGG_indexes = [];temp = [];

       %For relating a phyl dist between the modelled organism and the 
       %organisms on the BRENDA cell array it should search for a KEGG code
       %for each of these 
        for j=1:length(EC_indexes)

           %Assigns a KEGG index for those found on the KEGG struct
            orgs_index = find(strcmpi(orgs_cell(EC_indexes(j)),...
                                      TaxDistStruct.names),1);
            if ~isempty(orgs_index)
                KEGG_indexes = [KEGG_indexes; orgs_index];
                temp         = [temp;EC_indexes(j)];
           %For values related to organisms without KEGG code, then it will
           %look for KEGG code for the first organism with the same genus
            else
                k=1;
                while isempty(orgs_index) && k<length(TaxDistStruct.names)
                    str = TaxDistStruct.names{k};
                    org = orgs_cell{EC_indexes(j)};
                    if strcmpi(org(1:strfind(org,' ')-1),...
                            str(1:strfind(str,' ')-1))
                        orgs_index = k;KEGG_indexes = [KEGG_indexes;k];
                        temp = [temp;EC_indexes(j)];
                    end
                    k = k+1;
                end
            end
        end
       %Update the EC_indexes cell array
        EC_indexes = temp;
       %Looks for the taxonomically closest organism and saves the index 
       %of its appearences in the BRENDA cell
        if ~isempty(EC_indexes)
            distances  = num2cell(TaxDistStruct.distMat(org_index,:));
            distances  = distances(KEGG_indexes);
            EC_indexes = EC_indexes(cell2mat(distances) ==...
                                    min(cell2mat(distances)));                     
        end
    end
    
end
 
 
