function model = constrainModel(model, paramName, paramVal)
    for i = 1:length(paramName)
        if not(strcmp(paramName{i}, 'GROWTH'))
            model = setParam(model, 'lb', paramName{i}, paramVal(i));
            model = setParam(model, 'ub', paramName{i}, paramVal(i));
        end
    end
end

