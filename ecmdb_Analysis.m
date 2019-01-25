% % load ecmdb_ImportedDB_4-2018.mat
% % load Ecoli_iML1515_FVA_GeneDeletion.mat
% % 
% % % try to figure out if the KEGG IDs of metabolites in the EColi_iJO1366 are
% % % recorded in ecmdb
% % metsNum = length(mainModel.mets);
% % idStatus_KEGG = zeros(metsNum, 1);
% % idStatus_Pubchem = zeros(metsNum, 1);
% % idStatus_inchiKey = zeros(metsNum, 1);
% % 
% % for i = 1:metsNum
% %     if ~isempty(find(mainModel.EcoliKEGGIDs(i) == ecmdb.kegg_id_nums))
% %         idStatus_KEGG(i,1) = 1;
% %     end
% %     
% % end
% % 
% % for i = 1:metsNum
% %     if ~isempty(find(mainModel.pubChemIDs(i) == ecmdb.pubchem_id))
% %         idStatus_Pubchem(i,1) = 1;
% %     end 
% % end
% % 
% % for i = 1:metsNum
% %     if ~isempty(find(ismember(ecmdb.moldb_inchikey, mainModel.inchiKey{i}) == 1)) 
% %         idStatus_inchiKey(i,1) = 1;
% %     end
% % end
% 
% load ProductsDetails.mat
% load ProductsDetails_short.mat
% % prods = unique(prodIDsList);
% prods = (prodIDsList);
% load ecmdb_ImportedDB_4-2018.mat
% load Ecoli_iML1515.mat
% load SelectedOperators.mat
% keggID = zeros(length(prods),1);
% idStatus_Pubchem = zeros(length(prods),1);
% pubChemList = zeros(length(prods),1);
% for prIdx = 1:length(prods)
%     prIdx
%     if ~isempty(find(-1 * prods(prIdx) == ecmdb.pubchem_id))
%         idStatus_Pubchem(prIdx,1) = 1;
%         ecmdbKEGGId = find(-1 * prods(prIdx) == ecmdb.pubchem_id);
%         keggID(prIdx,1) = ecmdb.kegg_id_nums(ecmdbKEGGId(1));
%         pubChemList(prIdx,1) = -1 * prods(prIdx);
%     end
% end
% 
% keggID_unique = unique(keggID);
% pubChemList_unique = unique(pubChemList);
% keggIntersect = intersect(keggID_unique, mainModel.EcoliKEGGIDs);
% keggDiff = setdiff(keggID_unique, mainModel.EcoliKEGGIDs);
% keggIntersect(1) = [];
% 
% pubChemIntersect = intersect(pubChemList_unique, mainModel.pubChemIDs);
% prods = -1*prods;
% pubChemDiff = setdiff(prods, pubChemList_unique);
% pubChemIntersect(1) = [];
% 
% count = 1;
% details = [];
% for i = 1:length(keggDiff)
%     index = find(keggID == keggDiff(i));
%     for j = 1:length(index)
%         details = [details; substrateIDs(index(j)), prods(index(j)), keggDiff(i), enzymeOperatorIdx(index(j))];
%         count = count + 1;
%     end
% end
% prodMetFlag = zeros(length(mainModel.metNames),1);
% for i = 1:length(mainModel.metNames)
%     if ~isempty(find(mainModel.EcoliKEGGIDs(i) == keggIntersect)) ||...
%             ~isempty(find(mainModel.pubChemIDs(i) == pubChemIntersect))
%         prodMetFlag(i,1) = 1;
%         
%     end
% end
% 
% count = 1;
% pubChemDetails = [];
% detailsRxn = {};
% detailsEnzyme = {};
% for i = 1:length(pubChemDiff)
%     index = find(prods == pubChemDiff(i));
%     for j = 1:length(index)
%         pubChemDetails = [pubChemDetails; substrateIDs(index(j)), prods(index(j)), pubChemDiff(i), enzymeOperatorIdx(index(j))];
%         detailsRxn{count,1} = selectedOperators(enzymeOperatorIdx(index(j))).Reaction;
%         detailsEnzyme{count,1} = selectedOperators(enzymeOperatorIdx(index(j))).Enzyme;
%         count = count + 1;
%     end
% end

% 
% % load Ecoli_iML1515.mat
% % originalModel = mainModel;
% % 
% % load Ecoli_iML1515_FVA_GeneDeletion.mat
% % updatedModel = mainModel
% % 
% % for i = 1:length(updatedModel.metNames)
% %     index = find(ismember(originalModel.metNames, updatedModel.metNames{i}) == 1);
% %     if ~isempty(index) && (updatedModel.EcoliKEGGIDs(i) == 0)
% %         updatedModel.EcoliKEGGIDs(i) = originalModel.EcoliKEGGIDs(index(1));
% %     end
% % end
% 
% 
% % % find inchikey for all KEGG IDs in the model
% % load Ecoli_iML1515_FVA_GeneDeletion.mat
% % inchiKey = {};
% % 
% % for i = 1:length(mainModel.EcoliKEGGIDs)
% %     i
% %     keggID = mainModel.EcoliKEGGIDs(i)
% %     if keggID > 0
% %         keggIDStr = num2str(keggID);
% %         for j = 5-length(keggIDStr): -1 : 1
% %             keggIDStr = ['0', keggIDStr];
% %         end
% %         keggIDStr = ['C',keggIDStr];
% %         inchiKey{i,1} = char(py.getKEGGInfo.get_inchiKey(keggIDStr));
% %     else
% %         inchiKey{i,1} = '0';
% %     end
% % end
% 
% 
% % % find pubchemID for all mets in the model
% % % load Ecoli_iML1515_FVA_GeneDeletion.mat
% % load Ecoli_iML1515.mat
% % pubChemID = [];
% % 
% % for i = 1:length(mainModel.metNames)
% %     i
% %     metName = mainModel.metNames{i};
% %     pubID = str2num(char(py.OpenBabelFilesConversion.getCompoundsPubChemID(metName)));
% %     if isempty(pubID)
% %         pubChemID(i,1) = 0;
% %     else
% %         pubChemID(i,1) = pubID;
% %     end
% % end


load ecmdb_pubChemIDsLookup.mat
for i = 1:length(ecmdb.pubchem_idfromecmdb)
    if (ecmdb.pubchem_idfromecmdb(i) == ecmdb.PubChemIDfromInchiKeySara(i)) &&...
            
end
