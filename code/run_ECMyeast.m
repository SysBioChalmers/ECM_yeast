%function run_ECMyeast
%Update model files
cd model_curation
model = incorporateManualCuration;
clc
%Collect parameters and generate ECM_model
cd ../parameterize_model
ECM_model =  getParameterized_model(model);
%end