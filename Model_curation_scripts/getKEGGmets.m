function model = getKEGGmets(model,yeast)
% metsKEGG
%
% Function that gets a reduced model for yeast metabolism and maps the 
% all the metabolites into the concensus yeast8 model for the extraction of
% KEGG ID's for every metabolite.
%
% Ivan Domenzain 2018-07-11
%
metsList = cell(length(yeast.metNames),1);
metsKEGG = cell(length(model.metNames),1);
%Get list of original met strings
for i=1:length(yeast.metNames)
    metName = yeast.metNames{i};
    pos = strfind(metName,'[');
    if ~isempty(pos)
        metName = metName(1:pos(end)-2);
        metsList{i} =  metName;
    end
end

for i=1:length(model.metNames)
    metName = model.metNames{i};
    metName = nameExceptions(metName);
    
        
    index = find(strcmpi(metsList,metName),1);
    if ~isempty(index) && ~strcmpi(metName,'ubiquinone-6')
        KEGGID = yeast.metKEGGID{index};
        
        metsKEGG{i} = KEGGID;
    elseif strcmpi(metName,'HCO3')
    	metsKEGG{i} = 'C00288';
    elseif strcmpi(metName,'ubiquinone-6')
    	metsKEGG{i} = 'C00399';
    end

end
model.metKEGGID = metsKEGG;
end
function metName = nameExceptions(metName)
    if strcmpi(metName,'alpha-D-glucose')
        metName = 'D-glucose';
    elseif strcmpi(metName,'alpha-D-glucose 6-phosphate')
        metName = 'D-glucose 6-phosphate';
    elseif strcmpi(metName,'CO2')
        metName = 'carbon dioxide';  
    elseif strcmpi(metName,'O2')
        metName = 'oxygen';
    elseif strcmpi(metName,'D-glyceraldehyde 3-phosphate')
        metName = 'glyceraldehyde 3-phosphate';
    elseif strcmpi(metName,'D-ribose 5-phosphate')
        metName = 'ribose-5-phosphate';
    elseif strcmpi(metName,'glycerol-3-phosphate')
        metName = 'glycerol 3-phosphate';
    elseif strcmpi(metName,'glycerol-3-phosphate')
        metName = 'glycerol 3-phosphate';
    elseif strcmpi(metName,'nad+')
        metName = 'nad';
    elseif strcmpi(metName,'nadp+')
        metName = 'nadp(+)';    
    elseif strcmpi(metName,'3-phospho-D-glycerate')
        metName = '3-phosphonato-D-glycerate(3-)';
    elseif strcmpi(metName,'galactose')
        metName = 'D-galactose';
    elseif strcmpi(metName,'D-galactose 1-phosphate')
        metName = 'alpha-D-galactose 1-phosphate';
    elseif strcmpi(metName,'UDP-glucose')
        metName = 'UDP-D-glucose';
   elseif strcmpi(metName,'ubiquinol')
        metName = 'ubiquinol-6';
   elseif strcmpi(metName,'ubiquinone-9')
        metName = 'ubiquinone-6';
    end
end