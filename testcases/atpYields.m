addpath('sourceCode')
clf



tolerance = 0.01;

POratio = 0.9;


%Test changes stochiometry ATP synthase
% rxn = findIndex(model.rxns, 'ATP1');
% met1 = findIndex(model.mets, 'H[i]');
% met2 = findIndex(model.mets, 'H[m]');
% model.S(met1, rxn) = -3.3;
% model.S(met2, rxn) = 3.3;


%Load condtion
conditions = {  %'heyland2009';
                %'daran-lapujade2003';
                %'daran-lapujade2007';
                %'postma1989';
                %'hoek1999';
                'hoek1998';
                %'costenoble2011';
                'meyenburg1969';
            };
                    
load('../Models/model.mat');

%allow free exchange in absence of measurment
model = setParam(model,'ub','acOUT',0);
model = setParam(model,'ub','biomassOUT',1000);
model = setParam(model,'ub','co2OUT',1000);
model = setParam(model,'ub','ethOUT', 1000);
model = setParam(model,'ub','glyOUT', 1000);
model = setParam(model,'ub','PIOut', 1000);

model = setParam(model,'ub','acIN', 0);
model = setParam(model,'ub','o2IN',1000);
model = setParam(model,'ub','ethIN', 0);
model = setParam(model,'ub','PIIn', 1000);
model.b = [model.b model.b];

%Allow growth to be lower than measured
model = setParam(model, 'lb', 'GROWTH', 0);
model = setParam(model, 'ub', 'GROWTH', 1000);

%Set maintainance
model = setParam(model,'lb',{'ATPX'}, 0);  %mol/h maintainence            
model = setParam(model,'ub',{'ATPX'}, 1000);  %mol/h maintainence


valueObject = makeValueObjectWeight(0, 0, 0, 0, 0, 0, 0, 0, 1, 1);
model = makeBiomassEquation(model, valueObject);
model = setParam(model,'obj',{'GROWTH'}, 1);

[crap, exId] = getExchangeRxns(model);        
                
for i = 1:length(conditions)
    tmpModel = model;
    A = importdata(['conditions/' conditions{i} '.txt']);
    paramName = A.textdata;
    paramVal = A.data;
    
    growthRate = paramVal(findIndex(paramName, 'GROWTH'),:);
    if ismember('o2IN',paramName) == 1
        vO2 = paramVal(findIndex(paramName, 'o2IN'),:);
    else
        vO2 = paramVal(findIndex(paramName, 'co2OUT'),:);
        vO2 = vO2 -paramVal(findIndex(paramName, 'ethOUT'),:);
        vO2 = vO2 -paramVal(findIndex(paramName, 'acOUT'),:);
    end
    aerobicATP = vO2 .* (2 * POratio + 2/3);
    
    if ismember('ethOUT',paramName) == 1
        anaerobicATP = paramVal(findIndex(paramName, 'ethOUT'),:);
        if ismember('acOUT',paramName) == 1
            anaerobicATP = anaerobicATP + paramVal(findIndex(paramName, 'acOUT'),:);
            anaerobicATP = anaerobicATP - paramVal(findIndex(paramName, 'glyOUT'),:);
        end
    else
        deltaCO2 = paramVal(findIndex(paramName, 'co2OUT'),:)-vO2;
        anaerobicATP = deltaCO2;
    end
    totalATP = aerobicATP + anaerobicATP;
    
    subplot(1,2,1);
    title(['Calculated, P/O = ' num2str(POratio)])
    hold all
    plot(growthRate, totalATP, 'o-');
    ylim([0 50])
    xlabel('µ')
    ylabel('mmol ATP/gdw/h')
    
    subplot(1,2,2);
    title('FBA')
    hold all
    modelATP = zeros(length(growthRate),1);
    for j = 1:length(modelATP)
         tmpModel = constrainModel(tmpModel, paramName, paramVal(:,j), tolerance);
         tmpModel = setParam(tmpModel,'lb','glcIN',0);
         %tmpModel = setParam(tmpModel,'ub','co2OUT',1000);
         tmpModel = setParam(tmpModel,'lb','o2IN',0);
         tmpModel = setParam(tmpModel,'ub','glyOUT', 1000);
         tmpModel = setParam(tmpModel,'lb','GROWTH',0);
         tmpModel = setParam(tmpModel,'ub','GROWTH',1000);
         output = solveLin(tmpModel, true);
         modelATP(j) = -output.f;
    end
    plot(growthRate, modelATP, 'o-');
    ylim([0 50])
    xlabel('µ')
    ylabel('mmol ATP/gdw/h')

%     subplot(2,2,3);
%     hold all
%     plot(modelATP, totalATP, 'o-');
% 
%     subplot(2,2,4);
%     hold all
%     plot((modelATP-maintainance)./growthRate', totalATP, 'o-');
end

legend(conditions, 'location', 'NE')






