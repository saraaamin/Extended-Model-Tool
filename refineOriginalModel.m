% This function use Flux Variability Analysis from CobraToolbox to refine
% the model by removing all reactions with min and max flux = 0 when
% optimizing biomass production

function newModel = refineOriginalModel(mainModel)
% Get min and max flux of all reactions in the model
[minFlux, maxFlux] = fluxVariability(mainModel);

% minFlux = mainModel.minFlux;
% maxFlux = mainModel.maxFlux;

newModel = mainModel;
rxnsNum = length(minFlux);

% remove reactions with no flux going through them 
for i = 1:rxnsNum
    if minFlux(i) == 0 && maxFlux(i) == 0
        newModel = removeRxns(newModel, mainModel.rxns{i});
    end
end

% Remove genes that are no longer mapped to any reactions 
newModel = removeUnusedGenes(newModel);
end