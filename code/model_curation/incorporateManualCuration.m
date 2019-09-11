function model = incorporateManualCuration
model = importExcelModel('../../models/model.xls');
save('../../models/model.mat', 'model')
exportModel(model,'../../models/model.xml',true);
end
