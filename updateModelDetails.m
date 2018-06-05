function mainModel = updateModelDetails(mainModel)
% % get metabolites KEGG IDs
metsNum = length(mainModel.mets);
EcoliKEGGIDs = zeros(metsNum,1);

for i = 1:metsNum
    i
    if ismember('?',mainModel.metNames{i})
        EcoliKEGGIDs(i,1) = 0;
    else
        metKEGGID = py.getKEGGInfo.get_metabolite_ID(mainModel.metNames{i});
        EcoliKEGGIDs(i,1) = metKEGGID;
    end
end

% update the model with mets KEGG IDs
mainModel.EcoliKEGGIDs = EcoliKEGGIDs;

% % get ECnumbers from KEGG by use gene names
% ECNumList = [];
% genesNum = length(mainModel.genes);
% for i = 1:genesNum
%     i
%     ecNum = char(py.getKEGGInfo.get_EC_num(mainModel.genes{i}));
%     ecs = strsplit(ecNum, ',');
%     ECNumList{i,1} = ecs;
% end
% 
% mainModel.ECNumList = ECNumList;

% count = 1;
% for i = 1:length(mainModel.ECNumList)
%     if ~strcmp(mainModel.ECNumList{i}, '')
%         for j = 1:length(mainModel.ECNumList{i})
%             ECNumList{count,1} = mainModel.ECNumList{i}{j};
%             count = count + 1;
%         end
%     end
% end
% ECNumListUnique = unique(ECNumList);
% mainModel.ECNumListUnique = ECNumListUnique;

end