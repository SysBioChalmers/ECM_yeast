function model = buildKEGGrxns(model)
[~,n] = size(model.S);

KEGGrxns = cell(n,1);
for i=1:n
    rev = model.rev(i);
    keggRxn = '';
    RXNstr  = '';
    if ~isempty(model.grRules{i})
        substrates = find(model.S(:,i)<0);
        products   = find(model.S(:,i)>0);
        for j=1:length(substrates)
            index = substrates(j);
            coeff = abs(model.S(index,i));
            metKEGG = model.metKEGGID{index};
            
            if coeff~=1
                RXNstr = [RXNstr num2str(coeff) ' ' metKEGG];
            else
                RXNstr = [RXNstr metKEGG];
            end
            
            if j<length(substrates)
               RXNstr = [RXNstr ' + ']; 
            else
                RXNstr = [RXNstr ' '];
            end
        end
        
        if rev
            RXNstr = [RXNstr '<=> '];
        else
            RXNstr = [RXNstr '=> '];
        end
        
        for j=1:length(products)
            index   = products(j);
            coeff   = model.S(index,i);
            metKEGG = model.metKEGGID{index};
            
            if coeff~=1
                RXNstr = [RXNstr num2str(coeff) ' ' metKEGG];
            else
                RXNstr = [RXNstr metKEGG];
            end

            if j<length(products)
               RXNstr = [RXNstr ' + ']; 
            else
                RXNstr = [RXNstr ' '];
            end
        end
        
        disp(RXNstr)
    end
    KEGGrxns{i} = RXNstr;
end
model.KEGGrxns = KEGGrxns;
end