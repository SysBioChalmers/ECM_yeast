addpath('sourceCode')
load('../Model/model.mat');

%Load condtion
conditionName = 'heyland2009';
A = importdata(['conditions/' conditionName '.txt']);
paramName = A.textdata;
paramVal = A.data;

%allow free exchange in absence of measurment
model = setParam(model,'ub','o2IN',1000);
model = setParam(model,'ub','ethOUT', 1000);

%Set maintainance
model = setParam(model,'lb',{'ATPX'}, 1);  %0.7 mol/h maintainence

%Setup biomass equation
valueObject = makeValueObjectWeight(0.40, 0.12, 0.006, 0.025, 0.005, 0.01, 0, 0.4, 55, 1);
model = makeBiomassEquation(model, valueObject);

%Allow growth to be lower than measured
model = setParam(model, 'lb', 'GROWTH', 0);
model = setParam(model, 'ub', 'GROWTH', 1000);

%Maximize growth
model = setParam(model,'obj',{'GROWTH'}, 1);



results = zeros(length(model.rxns), size(paramVal,2));
for i = 1:size(results,2)
    model = constrainModel(model, paramName, paramVal(:,i));
    output = solveLinMin(model);
    if output.f == 0
        %if infeasible relax glucose uptake rate
        model = setParam(model, 'lb', 'glcIN', 0);
        output = solveLinMin(model);
    end
    results(:,i) = output.x;
end

subplot(1,2,1)
hold all
plot(paramVal(1,:), results(findIndex(model.rxns, 'GROWTH'),:), 'o')
plot([0 0.5], [0 0.5], 'k-')
xlim([0 0.5])
xlabel('measured growth')
ylabel('predicted growth')
axis square

subplot(1,2,2)
hold all
plot(paramVal(2,:), results(findIndex(model.rxns, 'glcIN'),:), 'o')
plot([0 25], [0 25], 'k-')
xlim([0 25])

xlabel('measured glucose')
ylabel('fitted glucose')
axis square

fileID = fopen(['fluxes/' conditionName '.txt'],'w');

fprintf(fileID, 'Rxn\tFlux\n');
for i = 1:length(resX)
    fprintf(fileID, '%s', model.rxns{i}); %Rxn Name
    for j = 1:size(results,2)
        fprintf(fileID, '\t%2.2f', results(i,j)); %Rxn flux
    end
    fprintf(fileID,'\n');
end
fclose(fileID);