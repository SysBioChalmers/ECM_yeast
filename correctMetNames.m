function model = correctMetNames(model)
    for l=1:length(model.metNames)
        if strcmpi(model.metNames(l),'Matrix protons')
            model.metNames(l) = {'h+'};
        elseif strcmpi(model.metNames(l),'free protons')
            model.metNames(l) = {'h+'};
        elseif strcmpi(model.metNames(l),'mitocondrialATP')
            model.metNames(l) = {'ATP'};
        elseif strcmpi(model.metNames(l),'mitocondrialADP')
            model.metNames(l) = {'ADP'};
        elseif strcmpi(model.metNames(l),'mitocondrial phosphate')
            model.metNames(l) = {'phosphate'};
        elseif strcmpi(model.metNames(l),'mitocondrialpyruvate')
            model.metNames(l) = {'pyruvate'};
        elseif strcmpi(model.metNames(l),'nad(+)')
            model.metNames(l) = {'nad+'};
        elseif strcmpi(model.metNames(l),'nadp(+)')
            model.metNames(l) = {'nadp+'};
        elseif strcmpi(model.metNames(l),'galactos')
            model.metNames(l) = {'galactose'};
        end
    end
end