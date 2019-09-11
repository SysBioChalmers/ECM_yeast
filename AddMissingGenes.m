%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function model = AddMissingGenes(model)
%
% Ivan Domenzain.   Last edited: 2018-02-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = AddMissingGenes(model)
    %Galactose metabolism genes that were lacking in the reduced model
    genes2add.genes   = {'YBR020W' 'YBR019' 'YBR018' 'YMR105C' 'YKL127'};
    genes2add.grRules = {'YBR020W' 'YBR019' 'YBR018' '(YMR105C or YKL127)'};
    genes2add.rxns    = [14 15 16 17];

    model = addGenes(model,genes2add.genes);
    for i=1:length(genes2add.grRules)
        model.grRules(genes2add.rxns(i)) =  genes2add.grRules(i);
    end
end


