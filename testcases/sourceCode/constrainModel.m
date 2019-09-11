function model = constrainModel(model, paramName, paramVal)
tolerance = 10^-3;

factor = 1-tolerance;
    for i = 1:length(paramName)
        model = setParam(model, 'lb', paramName{i}, paramVal(i)*factor);
        model = setParam(model, 'ub', paramName{i}, paramVal(i));
    end
end

