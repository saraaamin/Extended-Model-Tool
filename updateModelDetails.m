function mainModel = updateModelDetails(mainModel)
% % get metabolites KEGG IDs
metsNum = length(mainModel.mets);
EcoliKEGGIDs = zeros(metsNum,1);

for i = 1:metsNum
    if ismember('?',mainModel.metNames{i})
        EcoliKEGGIDs(i,1) = 0;
    else
        metKEGGID = py.getKEGGInfo.get_metabolite_ID(mainModel.metNames{i});
        EcoliKEGGIDs(i,1) = metKEGGID;
    end
end

% update the model with mets KEGG IDs
mainModel.EcoliKEGGIDs = EcoliKEGGIDs;

end