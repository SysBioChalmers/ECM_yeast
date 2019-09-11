addpath('sourceCode')
load('../Models/model.mat');
clf

%Load condtion
%conditionName = 'heyland2009';
%conditionName = 'costenoble2011';
%conditionName = 'daran-lapujade2003';
%conditionName = 'daran-lapujade2007';
conditionName = 'postma1989';
%conditionName = 'meyenburg1969';
%conditionName = 'hoek1999';
%conditionName = 'hoek1998';


A = importdata(['conditions/' conditionName '.txt']);
paramName = A.textdata;
paramVal = A.data;

%allow free exchange in absence of measurment
model = setParam(model,'ub','acOUT',0);
model = setParam(model,'ub','biomassOUT',1000);
model = setParam(model,'ub','co2OUT',1000);
model = setParam(model,'ub','ethOUT', 1000);
model = setParam(model,'ub','glyOUT', 1000);
model = setParam(model,'ub','glyOUT', 1000);
model = setParam(model,'ub','PIOut', 1000);

model = setParam(model,'ub','acIN', 0);
model = setParam(model,'ub','glcIN', 1000);
model = setParam(model,'ub','o2IN',1000);
model = setParam(model,'ub','ethIN', 0);
model = setParam(model,'ub','PIIn', 1000);

%Set maintainance
model = setParam(model,'lb',{'ATPX'}, 0.5);  %0.7 mol/h maintainence


%Setup biomass equation
valueObject = makeValueObjectWeight(0.40, 0.12, 0.006, 0.025, 0.005, 0.01, 0, 0.4, 30, 1);
%valueObject = makeValueObjectWeight(0, 0, 0, 0, 0, 0, 0, 0, 55, 1);
model = makeBiomassEquation(model, valueObject);



%Maximize growth
model = setParam(model,'obj',{'GROWTH'}, 1);



results = zeros(length(model.rxns), size(paramVal,2));
for i = 1:size(results,2)
    model = constrainModel(model, paramName, paramVal(:,i));

    %Allow growth to be lower than measured
    model = setParam(model, 'lb', 'GROWTH', 0);
    model = setParam(model, 'ub', 'GROWTH', 1000);
    
    output = solveLinMin(model, true);

    if output.f == 0
        fprintf('relaxing oxygen (%2.0f)\n', i)
        %if infeasible relax oxygen uptake growth
        model = constrainModel(model, paramName, paramVal(:,i));
        model = setParam(model,'lb','o2IN',0);
        model = setParam(model,'ub','o2IN',1000);

        output = solveLinMin(model, true);
        
        %if still infeasible relax glucose uptake rate
        if output.f == 0
            fprintf('relaxing glucose (%2.0f)\n', i)
            model = setParam(model, 'lb', 'glcIN', 0);
            output = solveLinMin(model, true);
        end
    end
    results(:,i) = output.x;
end


xlabels = {'predicted growth'
           'fitted glucose'
           'fitted O2'
           'fitted CO2'
           };
       
ylabels = {'measured growth'
          'measured glucose'
          'measured O2'
          'measured CO2'
    };
       
conditionRxns = {'GROWTH'
                 'glcIN'
                 'o2IN'
                 'co2OUT'
            };
       
for i = 1:4       
subplot(2,2,i)
hold all
curRxn = conditionRxns{i};

    if ismember(curRxn,paramName) == 1
        xVals = results(findIndex(model.rxns, curRxn),:);
        yVals = paramVal(findIndex(paramName,curRxn),:);

        plot(xVals, yVals, 'o')
        vmax = max(max(xVals), max(xVals));
        plot([0 vmax*1.1], [0 vmax*1.1], 'k-')
        xlim([0 vmax*1.1])
        ylim([0 vmax*1.1])
        xlabel(xlabels{i})
        ylabel(ylabels{i})
    end
axis square    
end


fileID = fopen(['fluxes/' conditionName '.txt'],'w');

fprintf(fileID, 'Rxn\tFlux\n');
for i = 1:length(model.rxns)
    fprintf(fileID, '%s', model.rxns{i}); %Rxn Name
    for j = 1:size(results,2)
        fprintf(fileID, '\t%2.2f', results(i,j)); %Rxn flux
    end
    fprintf(fileID,'\n');
end
fclose(fileID);