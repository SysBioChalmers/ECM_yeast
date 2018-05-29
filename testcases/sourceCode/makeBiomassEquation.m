function model = makeBiomassEquation(model, valueObject)
   biomassLocations.protein = findIndex(model.mets, 'growthProtein[c]');
   biomassLocations.rna = findIndex(model.mets, 'growthRNA[c]');
   biomassLocations.dna = findIndex(model.mets, 'growthDNA[c]');
   biomassLocations.lipid = findIndex(model.mets, 'growthLipid[c]');
   biomassLocations.glycogen = findIndex(model.mets, 'growthGlycogen[c]');
   biomassLocations.trehalose = findIndex(model.mets, 'growthTrehalose[c]');
   biomassLocations.mannan = findIndex(model.mets, 'growthMannan[c]');
   biomassLocations.glucan = findIndex(model.mets, 'growthGlucan[c]');
   biomassLocations.maintain = findIndex(model.mets, 'growthMaintainance[c]');
   biomassLocations.biomass = findIndex(model.mets, 'BIOMASS[c]');
   
   biomassEquation = findIndex(model.rxns, 'GROWTH');
   growthCol = makeGrowthCol(size(model.S, 1), biomassLocations, valueObject);
   model.S(:,biomassEquation) = growthCol;
end

