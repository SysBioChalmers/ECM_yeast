function model = mapRxnInfo(model,yeast)
% mapRxnInfo
%
% Function that maps every reaction in a yeast metabolic model and maps it
% into the concensus yeast8 model for extracting all the rxn-related
% curated information.
%
% Ivan Domenzain 2018-07-10
%
KEGGIDs = cell(length(model.rxns),1);
grRules = cell(length(model.rxns),1);
for i=1:length(model.rxnNames)
    rxn    = model.rxnNames{i};
    grRule = model.grRules{i};
    index  = find(strcmpi(yeast.rxnNames,rxn),1);
    KEGGID = '';
    if ~isempty(index)
        KEGGID  = yeast.rxnKEGGID{index};
        grRule  = '';
        %originalStr = yeast.rules{index};
        %grRule = createGrRule(originalStr,yeast,model,grRule);
    end
    %grRules{i} = grRule;
    KEGGIDs{i} = KEGGID;
end
%model.grRules = grRules;
[grRules,rxnGeneMat] = standardizeGrRules(model,false);
model.grRules = grRules;
model.rxnGeneMat = rxnGeneMat;
model.rxnKEGGID = KEGGIDs;
end

function grRule = createGrRule(originalStr,yeast,model,grRule)
if ~isempty(originalStr)
    %split the grRule into isoenzymes
    originalStr = strsplit(originalStr,'|');
    nIso = length(originalStr);
    for j=1:length(originalStr)
        isoenzyme = originalStr{j};
        %Remove all unused characters
        isoenzyme = strrep(isoenzyme,'x(','');
        isoenzyme = strrep(isoenzyme,'(','');
        isoenzyme = strrep(isoenzyme,')','');
        isoenzyme = strsplit(isoenzyme,' ');
        nUnits = length(isoenzyme);
        % get the isoenzyme genes
        for k=1:nUnits
            unit = isoenzyme{k};
            unit = strtrim(unit);
            %get the gene IDs for each numeric string
            if ~isempty(unit) && ~strcmpi(unit,'&')
                unit = str2num(unit);
                gene = yeast.genes{unit};
                if ~ismember(gene,model.genes)
                    model.genes = [model.genes;gene];
                end
                grRule = [grRule gene];
                % concatenate the different subunits with and's
            elseif strcmpi(unit,'&')
                grRule = [grRule 'and'];
            end
            %Add spaces if the isoenzyme hasn't finished yet
            if k<nUnits && nUnits>1
                grRule = [grRule ' '];
            end
        end
        if j<nIso && nIso>1
            grRule = [grRule ' or '];
        end
    end
end
end