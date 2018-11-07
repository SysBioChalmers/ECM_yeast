function model = correctMetNames(model)
current = pwd;
cd Databases


fileID = fopen('chem_prop.txt','r');
data = textscan(fileID,'%s %s %s %s %s %s %s %s %s','delimiter','\t');
fclose(fileID);
metaNetNames = data{2}(367:end);%(367:end);
InchiIds     = data{6}(367:end);
metNames     = model.metNames;%cell(length(model.mets,1));
i = 0;
    for l=1:length(model.metNames)
       if l~=4 
        metName = model.metNames{l};
        inchi   = model.inchis{l};
        inchi   = ['InChI=' inchi];
        inchi   = strrep(inchi,'InChI=1/','InChI=1S/');
        for j=1:length(InchiIds)
            datum = InchiIds{j};
            if strcmpi(datum,inchi)
                metNames{l} = metaNetNames{j};
                i = i+1;
            end
        end
%         disp(inchi)
%         index   = find(strcmpi(InchiIds,inchi),1);
%         if ~isempty(index)
%             i=i+1;
%             metNames{index} = metaNetNames{index};
%         end
        %metNames{l} = manualAssignment(metName);
        %metNames{l} = metNames{l};%(1);
        
        if strcmpi(metNames(l),'galactos')
            metNames{l} = 'D-galactose';
        elseif strcmpi(metNames(l),'1,3-bis-phosphoglycerate')
            metNames{l} = '1,3-bisphosphoglycerate';
        end
        
        metNames{l} = char(metNames{l});%metaNetNames{j};
        disp(['Ready with metabolite: ' metNames{l}  ' ' num2str(i)])
       end
    end
    model.metNames = metNames;
cd (current)
end

function metName = manualAssignment(metName)
        if strcmpi(metName,'Matrix protons')
            metName = {'h+'};
        elseif strcmpi(metName,'free protons')
            metName = {'h+'};
        elseif strcmpi(metName,'mitocondrialATP')
            metName = {'ATP'};
        elseif strcmpi(metName,'mitocondrialADP')
            metName = {'ADP'};
        elseif strcmpi(metName,'mitocondrial phosphate')
            metName = {'phosphate'};
        elseif strcmpi(metName,'mitocondrialpyruvate')
            metName = {'pyruvate'};
        elseif strcmpi(metName,'nad(+)')
            metName = {'nad+'};
        elseif strcmpi(metName,'nadp(+)')
            metName = {'nadp+'};
        elseif strcmpi(metName,'galactos')
            metName = {'galactose'};
        elseif strcmpi(metName,'1,3-bis-phosphoglycerate')
            metName = '1,3-bisphosphoglycerate';
        end
end