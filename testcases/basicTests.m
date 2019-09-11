addpath('sourceCode')
load('../Models/model.mat');

%allow free exchange in absence of measurment
model = setParam(model,'ub','acOUT',0);
model = setParam(model,'ub','biomassOUT',1000);
model = setParam(model,'ub','co2OUT',1000);
model = setParam(model,'ub','ethOUT', 1000);
model = setParam(model,'ub','glyOUT', 1000);
model = setParam(model,'ub','glyOUT', 1000);
model = setParam(model,'ub','PIOut', 1000);

model = setParam(model,'ub','acIN', 0);
model = setParam(model,'ub','o2IN',1000);
model = setParam(model,'ub','ethIN', 0);
model = setParam(model,'ub','PIIn', 1000);
model = setParam(model,'ub','HCO3IN', 1000);

%Set maintainance

model.b = [model.b model.b];

%Setup biomass equation
%valueObject = makeValueObjectWeight(0.40, 0.12, 0.006, 0.025, 0.005, 0.01, 0, 0.4, 55, 1);


%Allow growth to be lower than measured
model = setParam(model, 'lb', 'GROWTH', 0);
model = setParam(model, 'ub', 'GROWTH', 1000);

%Maximize ATP
model = setParam(model,'ub','glcIN', 1);
model = setParam(model,'lb',{'ATPX'}, 0);  %0.7 mol/h maintainence
valueObject = makeValueObjectWeight(0, 0, 0, 0, 0, 0, 0, 0, 1, 1);
model = makeBiomassEquation(model, valueObject);
model = setParam(model,'obj',{'GROWTH'}, 1);
output = solveLinMin(model, true);
-output.f/1

%Maximize biomass
model = setParam(model,'ub','glcIN', 1);
model = setParam(model,'lb',{'ATPX'}, 0.5);  %0.7 mol/h maintainence
valueObject = makeValueObjectWeight(0.40, 0.12, 0.006, 0.025, 0.005, 0.01, 0, 0.4, 55, 1);
model = makeBiomassEquation(model, valueObject);
model = setParam(model,'obj',{'GROWTH'}, 1);
output = solveLinMin(model, true);
-output.f/0.180

%Maximize biomass on ethanol
model = setParam(model,'ub','glcIN', 0);
model = setParam(model,'ub','ethIN', 1);
model = setParam(model,'lb',{'ATPX'}, 0.5);  %0.7 mol/h maintainence
valueObject = makeValueObjectWeight(0.40, 0.12, 0.006, 0.025, 0.005, 0.01, 0, 0.4, 55, 1);
model = makeBiomassEquation(model, valueObject);
model = setParam(model,'obj',{'GROWTH'}, 1);
output = solveLinMin(model, true);
-output.f/0.046

%Maximize biomass on acetate
model = setParam(model,'ub','ethIN', 0);
model = setParam(model,'ub','acIN', 1);
model = setParam(model,'lb',{'ATPX'}, 0.5);  %0.7 mol/h maintainence
valueObject = makeValueObjectWeight(0.40, 0.12, 0.006, 0.025, 0.005, 0.01, 0, 0.4, 55, 1);
model = makeBiomassEquation(model, valueObject);
model = setParam(model,'obj',{'GROWTH'}, 1);
output = solveLinMin(model, true);
-output.f/0.060
