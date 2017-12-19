%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Kcats,] = ParametersCoverage(model_data)
%
% Kinetic Parameters coverage estimation (Kcats and Km values for substrates
% and products)for all of the enzymatic reactions  with an associated EC 
% number in model_data.model. 
%
% Ivan Domenzain.   Last edited: 2017-12-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EC_numbers    = model_data.EC_numbers;
model         = model_data.model;
[nR, x]       = size(EC_numbers);
Totalcounter  = 0;
Kcatcounter  = 0;
KScounter     = 0;
KPcounter     = 0;

for i=1:nR
    subsIndx      = find(model.S(:,i)<0);
    prodsIndx     = find(model.S(:,i)>0);
    substrates    = model.metNames(subsIndx);
    products      = model.metNames(prodsIndx);
    nonEmpty      = find(~cellfun(@isempty, EC_numbers(i,:)));
    rxnECs        = unique(EC_numbers(i,nonEmpty));
    Kcatcounter   = Kcatcounter + length(rxnECs);
    KScounter     = KScounter + length(rxnECs)*(length(substrates));
    KPcounter     = KPcounter + length(rxnECs)*(length(products));
    Totalcounter  = Totalcounter + length(rxnECs) + ...
                    length(rxnECs)*(length(substrates) + length(products));
                
                
    % Next step! to read data from 
                      
end

