% This function predicts derivatives of enzyme promiscuity 
function extendModel(mainModel)
global substrateIDs
global prodIDsList
global formulaList

originalModel = mainModel;
originalEColiKeggID = mainModel.EcoliKEGGIDs;

rxnIDList = [];
productsFolder = '.\productsMol\';
savedProductsFolder = '.\savedProductsMol\';
delete([savedProductsFolder, '\*.mol']);

uniqueEColiKeggID = unique(originalEColiKeggID);
cids = [];
cidNames = [];
cidCount = 0;
unknownProdIDs = [];
unknownFormulasList = [];
substrateIDs = [];
prodIDsList = [];
formulaList = [];
enzymeOperatorIdx = [];
prodEnzymesList = [];
stepsFBAResults = [];
allModel_Combined=[];
originalRxn_couplingList = [];
operatorsPerStepList = [];
newProdID = -101.01;

% set of operators associated with the enzymes of the model
operatorsFileName = 'SelectedOperators';
load(operatorsFileName)

% metabolites with high concentration
load cofactors.mat
compoundsToApplyOpsOn = [3; 6; 9; 19; 22; 24; 25; 26; 29; 31; 37; 41;...
    42; 43; 47; 49; 51; 53; 54; 58; 59; 62; 64; 65; 73; 74; 78; 79; 82;...
    85; 91; 92; 100; 106; 108; 111; 118; 119; 122; 127; 135; 141; 147;...
    148; 149; 152; 155; 156; 158; 166; 183; 188; 196; 197; 199; 227; 234;...
    236; 242; 257; 258; 279; 288; 299; 318; 327; 332; 337; 345; 352; 354;...
    380; 387; 407; 417; 437; 438; 469; 475; 493; 559; 620; 631; 666; 672;...
    860; 1236; 1602; 2504; 2631; 3175; 3722; 3736; 3794; 4256; 4376; 4411;...
    4677; 4823; 4874; 5382; 5754; 5809; 6022; 17556; 17569];

% compound IDs that don't have matches in KEGG or compounds consisting of
% one atom with no bounds. Those compounds are skipped since no products 
% can be generated from them
% the last set of IDs are carrier proteins that should be excluded
% since they contain S and R groups which can't be balanced as they are
% not detailed.
    excludedIDs = [-1; 0; 1; 10; 14; 23; 34; 38; 70; 76; 80; 84; 87; 175; 238; 282; 283; 291; 305; 698; 703; 787; 824; 1330; 1342; 1413; 1528; 1635;...
        1636; 1637; 1638; 1639; 1640; 1641; 1642; 1643; 1644; 1645; 1646; 1647; 1648; 1649; 1650; 1651; 1652; 1653; 1834; 2386; 2745; 2869;...
        5737; 6710; 7292; 14818; 14819; 15233; 19610;...
        173; 4180; 4619; 4620; 4688; 5223; 5274];

% for each compound in the model, apply the operators to generate products
% associated to the operators
prodCounter = 1;

for compoundIdx = 1:length(compoundsToApplyOpsOn)
    compound = compoundsToApplyOpsOn(compoundIdx);

    if ~isempty(find(excludedIDs == compound))
        continue;
    end
    
    % get the string of a compound ID to look up the KCF file on KEGG
    digitsNum = floor(log10 (compound))+1;
    f = '';
    if digitsNum~=5
        f = '0';
        for i = 1:5-1-digitsNum
            f = strcat(f,'0');
        end
    end
    
    compoundstr = strcat('C',f,num2str(compound));
    
    % generate a list of products when applying the operators to the
    % current compound 
    [inputList, inputStructure, inputListNumbering, products, operatorIdx] = GenerateProducts(compoundstr,operatorsFileName);
    
    % extend the model if only there are products predicted by Proximal
    if ~isempty(products)
        [prodEnzymes, operatorIdx] = GenerateMolFiles(inputStructure, inputListNumbering, products, operatorIdx);
        
        productsFolderFiles = dir(productsFolder);
        
        for prodIdx = 3:length(productsFolderFiles)
            prodIdx    
            prodFormula = [];
            currentProdFile = productsFolderFiles(prodIdx).name;
            fileIdx = str2double(productsFolderFiles(prodIdx).name(9:end-4));
            prodKCFFile = generate_kcfFile([productsFolder , currentProdFile]);
            prodCmpID = compareProdKCFToKeggKCF(prodKCFFile);
%             prodCmpID =  [];
            
            if prodCmpID == compound
                continue;
            end
                
            operatorsPathwayDetailsIdx = operatorIdx(fileIdx);
            
            if ~isempty(prodCmpID)
                % keeping track of compounds and products without
                % adding them to the model
                if length(prodCmpID) == 1
                    switchCase = 1;
                    substrateIDs = [substrateIDs; compound];
                    prodIDsList = [prodIDsList; prodCmpID];
                    enzymeOperatorIdx = [enzymeOperatorIdx; operatorsPathwayDetailsIdx];
                    prodEnzymesList{prodCounter,1} = prodEnzymes{fileIdx};
                    prodCounter = prodCounter + 1;
                else
                    switchCase = 2;
                    for index = 1:length(prodCmpID)
                        substrateIDs = [substrateIDs; compound];
                        prodIDsList = [prodIDsList; prodCmpID(index)];
                        enzymeOperatorIdx = [enzymeOperatorIdx; operatorsPathwayDetailsIdx];
                        prodEnzymesList{prodCounter,1} = prodEnzymes{fileIdx};
                        prodCounter = prodCounter + 1;
                    end
                end
            else
                % retrieve compound details using OpenBabel toolbox
                switchCase = 1;
                productOBDetails = py.OpenBabelFilesConversion.main_function(currentProdFile);
                
                if ~isempty(productOBDetails) 
                    prodOBStr = char(py.str(productOBDetails));
                    separatorIdx = find(prodOBStr == ',');
                    
                    if ~isempty(separatorIdx)
                        pubChemID = str2double(prodOBStr(1:separatorIdx-1));
                        
                        if ~isnan(pubChemID)
                            % setting PubChem IDs to -ive to differentiate between them
                            % and IDs retreived by KEGG
                            prodCmpID = -1*pubChemID;
                            pubChemName = prodOBStr(separatorIdx+1:end);
                            cidCount = cidCount + 1;

                            if isempty(find(prodCmpID == cids))
                                cids(cidCount,1) = prodCmpID;
                                cidNames{cidCount,1} = pubChemName;
                                prodFormula = py.OpenBabelFilesConversion.convert_mol_to_formula(currentProdFile);
                                prodFormula = char(py.str(prodFormula));
                                formulaList = [formulaList, prodFormula];
                                
                                % keeping track of compounds and products without
                                % adding them to the model
                                substrateIDs = [substrateIDs; compound];
                                prodIDsList = [prodIDsList; prodCmpID];
                                enzymeOperatorIdx = [enzymeOperatorIdx; operatorsPathwayDetailsIdx];
                                prodEnzymesList{prodCounter,1} = prodEnzymes{fileIdx};
                                prodCounter = prodCounter + 1;
                            end
                        else
                            % continue;
                            
                            % if a product is not found in PubChem nor
                            % KEGG, assign an ID for it
                            prodCmpID = newProdID;
                            prodFormula = py.OpenBabelFilesConversion.convert_mol_to_formula(currentProdFile);
                            prodFormula = char(py.str(prodFormula));
                            unknownProdIDs(end+1,1) = prodCmpID;
                            unknownFormulasList{end+1, 1} = prodFormula;
                            newProdID = newProdID - 1;
                            
                            % keeping track of compounds and products without
                            % adding them to the model
                            substrateIDs = [substrateIDs; compound];
                            prodIDsList = [prodIDsList; prodCmpID];
                            % saving the mol file in a different directory for
                            % later consideration
                            newNameProdFile = [currentProdFile(1:8), num2str(compound), num2str(prodCmpID), '.mol'];
                            sourceFolder = [productsFolder, currentProdFile];
                            destinationFolder = [savedProductsFolder, newNameProdFile];
                            movefile(sourceFolder, destinationFolder);

                        end
                    else
                        % continue;
                        
                        % if a product is not found in PubChem nor
                        % KEGG, assign an ID for it
                        prodCmpID = newProdID;
                        prodFormula = py.OpenBabelFilesConversion.convert_mol_to_formula(currentProdFile);
                        prodFormula = char(py.str(prodFormula));
                        unknownProdIDs(end+1,1) = prodCmpID;
                        unknownFormulasList{end+1, 1} = prodFormula;
                        newProdID = newProdID - 1;
                        
                        % keeping track of compounds and products without
                        % adding them to the model
                        substrateIDs = [substrateIDs; compound];
                        prodIDsList = [prodIDsList; prodCmpID];
                        % saving the mol file in a different directory for
                        % later consideration
                        newNameProdFile = [currentProdFile(1:end-4), '-', num2str(compound), num2str(prodCmpID), '.mol'];
                        sourceFolder = [productsFolder, currentProdFile];
                        destinationFolder = [savedProductsFolder, newNameProdFile];
                        movefile(sourceFolder, destinationFolder);

                    end
                   
                else
                    % continue;
                    
                    % if a product is not found in PubChem nor
                    % KEGG, assign an ID for it
                    prodCmpID = newProdID;
                    prodFormula = py.OpenBabelFilesConversion.convert_mol_to_formula(currentProdFile);
                    prodFormula = char(py.str(prodFormula));
                    unknownProdIDs(end+1,1) = prodCmpID;
                    unknownFormulasList{end+1, 1} = prodFormula;
                    newProdID = newProdID - 1;
                    % keeping track of compounds and products without
                    % adding them to the model
                    substrateIDs = [substrateIDs; compound];
                    prodIDsList = [prodIDsList; prodCmpID];
                    % saving the mol file in a different directory for
                    % later consideration
                    newNameProdFile = [currentProdFile(1:end-4), '-', num2str(compound), num2str(prodCmpID), '.mol'];
                    sourceFolder = [productsFolder, currentProdFile];
                    destinationFolder = [savedProductsFolder, newNameProdFile];
                    movefile(sourceFolder, destinationFolder);

                end
                
            end  
            
            % If the identified product is one of the cofactors, ignore it
%             if ~isempty(find(gluCoFactKEGGID == prodCmpID(1)))
%                 continue;
%             end
%             
%             if isempty(operatorsPathwayDetailsIdx)
%                 continue;
%             end
            
%             if switchCase == 1
%                 index = find(originalEColiKeggID == prodCmpID);
%                 [updatedRxnIDList, updatedModel, modelAugmentFlag] = updateModel_OneStep(originalModel, originalEColiKeggID, index, compound, prodCmpID, selectedOperators(operatorsPathwayDetailsIdx));
% 
%                 
%                 if modelAugmentFlag
%                     % couple the reaction so the flux would be in the same
%                     % direction as the original reaction flux
%                     [updatedModel_PosRxn, updatedEColiKeggID, mainCouplingRxn_PosRxn, posBoundVal] = setBoundsForAddedRxn(updatedModel, originalEColiKeggID, updatedRxnIDList(end), couplePercent, 1, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1}, selectedOperators(operatorsPathwayDetailsIdx));
%                     mutantPosFluxValue = optimizeCbModel(updatedModel_PosRxn, 'max', false, false);
% 
%                     % couple the reaction so the flux would be in the
%                     % opposite direction as the original reaction flux
% %                     [updatedModel_NegRxn, updatedEColiKeggID, mainCouplingRxn_NegRxn] = coupleRxnFlux(updatedModel, originalEColiKeggID, prodEnzymesList{prodIdx-2,1},...
% %                         pathwayFluxDetails, updatedRxnIDList(end), couplePercent, -1);
% 
% %                     [updatedModel_NegRxn, updatedEColiKeggID, mainCouplingRxn_NegRxn] = restrictFluxOfAddedReaction(compound, updatedModel, originalEColiKeggID, updatedRxnIDList(end), couplePercent, -1, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1});
%                     [updatedModel_NegRxn, updatedEColiKeggID, mainCouplingRxn_NegRxn, negBoundVal] = setBoundsForAddedRxn(updatedModel, originalEColiKeggID, updatedRxnIDList(end), couplePercent, -1, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1}, selectedOperators(operatorsPathwayDetailsIdx));
%                     mutantNegFluxValue = optimizeCbModel(updatedModel_NegRxn, 'max', false, false);
% 
% 
%                     % setting default for copuling to be in the same direction
%                     coupleDirection = 1;
%                     mutantFluxValue = mutantPosFluxValue; 
%                     boundVal = posBoundVal;
%                     
%                     % check if the default coupling direction will change                
%                     if mutantPosFluxValue.f < mutantNegFluxValue.f && mutantPosFluxValue.f > 0.0001
%                         coupleDirection = 1;
%                         mutantFluxValue = mutantPosFluxValue; 
%                         boundVal = posBoundVal;
%                     elseif mutantNegFluxValue.f < mutantPosFluxValue.f && mutantNegFluxValue.f > 0.0001
%                         coupleDirection = -1;
%                         mutantFluxValue = mutantNegFluxValue;
%                         boundVal = negBoundVal;
%                     end
%                     
%                     modelAugmentFlag = false;
% 
%                     if mutantFluxValue.f > 0.0001
%                         index = find(EColiKeggID == prodCmpID);
%                         [rxnIDList, modelAugmentFlag] = updateMainModel(index, compound, prodCmpID, rxnIDList, selectedOperators(operatorsPathwayDetailsIdx));
%                     end
% 
%                     if modelAugmentFlag
%                         stepsFBAResultsPos_OneStep = [stepsFBAResultsPos_OneStep; mutantPosFluxValue];
%                         allModelPos_OneStep = [allModelPos_OneStep; updatedModel_PosRxn];
% 
%                         stepsFBAResultsNeg_OneStep = [stepsFBAResultsNeg_OneStep; mutantNegFluxValue];
%                         allModelNeg_OneStep = [allModelNeg_OneStep; updatedModel_NegRxn];
% 
% %                         [mainModel, EColiKeggID, mainCouplingRx] = restrictFluxOfAddedReaction(compound, mainModel, EColiKeggID, rxnIDList(end), couplePercent, coupleDirection, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1});
%                         [mainModel, EColiKeggID, mainCouplingRx] = setBoundsForAddedRxn_FinalModel(mainModel, EColiKeggID, rxnIDList(end), coupleDirection, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1}, boundVal, selectedOperators(operatorsPathwayDetailsIdx));
%                     
% %                         [mainModel, EColiKeggID, mainCouplingRxn] = coupleRxnFlux(mainModel, EColiKeggID,...
% %                             prodEnzymesList{prodIdx-2,1}, pathwayFluxDetails, rxnIDList(end), couplePercent, coupleDirection); 
%     %                     mutantFluxValue = optimizeCbModel(mainModel, 'max', false, false);
%     %                     stepsFBAResults = [stepsFBAResults; mutantFluxValue];
%     %                     originalRxn_couplingList = [originalRxn_couplingList; mainCouplingRxn];
% 
%                         selectedOperators(operatorsPathwayDetailsIdx).moleculePosition = products(fileIdx).ObtainedFromInput;
%                         selectedOperators(operatorsPathwayDetailsIdx).coupleDirection = coupleDirection;
%                         selectedOperators(operatorsPathwayDetailsIdx).OperatorIdx = operatorsPathwayDetailsIdx;
% 
%                         operatorsPerStepList = [operatorsPerStepList; selectedOperators(operatorsPathwayDetailsIdx)];
%                         allModel_Combined = [allModel_Combined; mainModel];

%                         % saving the mol file in a different directory for
%                         % later consideration
%                         newNameProdFile = [currentProdFile(1:end-4), '-', num2str(prodCmpID), '.mol'];
%                         sourceFolder = [productsFolder, currentProdFile];
%                         destinationFolder = [savedProductsFolder, newNameProdFile];
%                         movefile(sourceFolder, destinationFolder);
%                         
%                     end
% 
%                 end
%                 
%                 
%             elseif switchCase == 2
%                 currentCouplePercent = couplePercent/length(prodCmpID);
%                 
%                 for kfcProdIdx = 1:length(prodCmpID)
%                     if ~isempty(find(gluCoFactKEGGID == prodCmpID(kfcProdIdx),1))
%                         continue;
%                     end
% 
%                     index = find(originalEColiKeggID == prodCmpID(kfcProdIdx));
%                     [updatedRxnIDList, updatedModel, modelAugmentFlag] = updateModel_OneStep(originalModel, originalEColiKeggID, index, compound, prodCmpID(kfcProdIdx), selectedOperators(operatorsPathwayDetailsIdx));
%                    
%                     if modelAugmentFlag
%                         % couple the reaction so the flux would be in the same
%                         % direction as the original reaction flux
% %                         [updatedModel_PosRxn, updatedEColiKeggID, mainCouplingRxn_PosRxn] = coupleRxnFlux(updatedModel, originalEColiKeggID, prodEnzymesList{prodIdx-2,1},...
% %                                             pathwayFluxDetails, updatedRxnIDList(end), couplePercent, 1);
% %                         [updatedModel_PosRxn, updatedEColiKeggID, mainCouplingRxn_PosRxn] = restrictFluxOfAddedReaction(compound, updatedModel, originalEColiKeggID, updatedRxnIDList(end), couplePercent, 1, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1});
%                         [updatedModel_PosRxn, updatedEColiKeggID, mainCouplingRxn_PosRxn, posBoundVal] = setBoundsForAddedRxn(updatedModel, originalEColiKeggID, updatedRxnIDList(end), couplePercent, 1, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1}, selectedOperators(operatorsPathwayDetailsIdx));
%                         mutantPosFluxValue = optimizeCbModel(updatedModel_PosRxn, 'max', false, false);
% 
%                         % couple the reaction so the flux would be in the
%                         % opposite direction as the original reaction flux
% %                         [updatedModel_NegRxn, updatedEColiKeggID, mainCouplingRxn_NegRxn] = coupleRxnFlux(updatedModel, originalEColiKeggID, prodEnzymesList{prodIdx-2,1},...
% %                                             pathwayFluxDetails, updatedRxnIDList(end), couplePercent, -1);
% 
% %                         [updatedModel_NegRxn, updatedEColiKeggID, mainCouplingRxn_NegRxn] = restrictFluxOfAddedReaction(compound, updatedModel, originalEColiKeggID, updatedRxnIDList(end), couplePercent, -1, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1});
%                         [updatedModel_NegRxn, updatedEColiKeggID, mainCouplingRxn_NegRxn, negBoundVal] = setBoundsForAddedRxn(updatedModel, originalEColiKeggID, updatedRxnIDList(end), couplePercent, -1, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1}, selectedOperators(operatorsPathwayDetailsIdx));
%                         mutantNegFluxValue = optimizeCbModel(updatedModel_NegRxn, 'max', false, false);
% 
% 
%                         % setting default for copuling to be in the same direction
%                         coupleDirection = 1;
%                         mutantFluxValue = mutantPosFluxValue; 
%                         boundVal = posBoundVal;
%                         
%                         if (mutantPosFluxValue.f <= 0.0001 || mutantNegFluxValue.f <= 0.0001)
%                             i;
%                         end
%                     
%                         % check if the default coupling direction will change                
%                         if mutantPosFluxValue.f < mutantNegFluxValue.f && mutantPosFluxValue.f > 0.0001
%                             coupleDirection = 1;
%                             mutantFluxValue = mutantPosFluxValue; 
%                             boundVal = posBoundVal;
%                         elseif mutantNegFluxValue.f < mutantPosFluxValue.f && mutantNegFluxValue.f > 0.0001
%                             coupleDirection = -1;
%                             mutantFluxValue = mutantNegFluxValue;
%                             boundVal = negBoundVal;
%                         end
%                         modelAugmentFlag = false;
%                     
%                         if mutantFluxValue.f > 0.0001
%                             index = find(EColiKeggID == prodCmpID(kfcProdIdx));
%                             [rxnIDList, modelAugmentFlag] = updateMainModel(index, compound, prodCmpID(kfcProdIdx), rxnIDList, selectedOperators(operatorsPathwayDetailsIdx));
%                         end
% 
%                         if modelAugmentFlag
%                             stepsFBAResultsPos_OneStep = [stepsFBAResultsPos_OneStep; mutantPosFluxValue];
%                             allModelPos_OneStep = [allModelPos_OneStep; updatedModel_PosRxn];
% 
%                             stepsFBAResultsNeg_OneStep = [stepsFBAResultsNeg_OneStep; mutantNegFluxValue];
%                             allModelNeg_OneStep = [allModelNeg_OneStep; updatedModel_NegRxn];
%                             
% %                             [mainModel, EColiKeggID, mainCouplingRx] = restrictFluxOfAddedReaction(compound, mainModel, EColiKeggID, rxnIDList(end), couplePercent, coupleDirection, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1});
%                             [mainModel, EColiKeggID, mainCouplingRx] = setBoundsForAddedRxn_FinalModel(mainModel, EColiKeggID, rxnIDList(end), coupleDirection, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1}, boundVal, selectedOperators(operatorsPathwayDetailsIdx));
%                             
% %                             [mainModel, EColiKeggID, mainCouplingRxn] = coupleRxnFlux(mainModel, EColiKeggID,...
% %                                 prodEnzymesList{prodIdx-2,1}, pathwayFluxDetails, rxnIDList(end), currentCouplePercent, coupleDirection); 
%     %                         mutantFluxValue = optimizeCbModel(mainModel, 'max', false, false);
%     %                         stepsFBAResults = [stepsFBAResults; mutantFluxValue];
%     %                         originalRxn_couplingList = [originalRxn_couplingList; mainCouplingRxn];
% 
%                             selectedOperators(operatorsPathwayDetailsIdx).moleculePosition = products(prodIdx-2).ObtainedFromInput;
%                             selectedOperators(operatorsPathwayDetailsIdx).coupleDirection = coupleDirection;
%                             selectedOperators(operatorsPathwayDetailsIdx).OperatorIdx = operatorsPathwayDetailsIdx;
% 
%                             operatorsPerStepList = [operatorsPerStepList; selectedOperators(operatorsPathwayDetailsIdx)];
%                             
%                             allModel_Combined = [allModel_Combined; mainModel];
% 
%                             % saving the mol file in a different directory for
%                             % later consideration
%                             if kfcProdIdx == 1
%                                 newNameProdFile = [currentProdFile(1:end-4), '-', num2str(prodCmpID(kfcProdIdx)), '.mol'];
%                                 sourceFolder = [productsFolder, currentProdFile];
%                                 destinationFolder = [savedProductsFolder, newNameProdFile];
%                                 movefile(sourceFolder, destinationFolder);
%                             end
% 
%                         end
% 
%                     end
%                     
%                      
%                 end
%             end
%             end
            
        end
    end
end

save ProductsDetails_withunknown.mat prodIDsList substrateIDs formulaList enzymeOperatorIdx prodEnzymesList unknownProdIDs unknownFormulasList
end
