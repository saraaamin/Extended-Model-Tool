% As an input later the user can define the model to be extended
% As a requirement:
%   - should provide the data files mapping compounds and reactions in KEGG

function run_ExtendMetabolicModel()
% initCobraToolbox()
global reactions
global reactionData
global indexReactions
global indexCompounds
global compounds

% Files containing KEGG compounds and reactions mapping
load reactionData_012716.mat
load reactions_012716.mat
load indexReactions_012716.mat

load compounds_012716.mat
load indexCompounds_012716.mat

% Loading the model and it's KEGG IDs
load Ecoli_iML1515.mat

% Refine model using FVA and gene deletion
% mainModel = refineOriginalModel(mainModel);

%Get model compounds details
% mainModel = updateModelDetails(mainModel);

% Get unique enzymes from the model
% enzymeList = getEnzymeList(mainModel);

% for user given models use the EC numbers associated with genes in the 
% model to get the operators mapped to them
% selectedOperators = getOperators(enzymeList);    

extendModel(mainModel)
end

function selectedOperators = getOperators(enzymesList)
global operators
global operatorsMainAtom
global operatorsMainAtomPos

% load the Operators Neda generated from KEGG
% load KEGGOperators_2016_08_14.mat
% load KEGGOperators_2016_08_14_update.mat
load operators.mat

enzymeNum = length(enzymesList);

% for each enzyme, get the list of enzymes catalyzing the reactions
selectedOperators = struct();
% each enzyme in the enzymesList, get the list of operators associated
% to it from the allOperators dataset
for enzymeIdx = 1:enzymeNum
    enzymeIdx
%     selectedOperators = getEnzymeOperators(enzymesList{enzymeIdx}, rxnKEGGIDList(enzymeIdx), selectedOperators);
    selectedOperators = getEnzymeOperators(enzymesList{enzymeIdx}, selectedOperators);
end

selectedOperators = getUniqueOperators(selectedOperators);

[selectedOperators, operatorsMainAtom, operatorsMainAtomPos] = sortRs(selectedOperators);
save SelectedOperators selectedOperators operatorsMainAtom operatorsMainAtomPos


end

% This is a function that returns the operators associated with a specific
% enzyme.
% input:
%   enzyme: EC number of the quired enzyme
%   selectedOperators: a list of current selected operators in case the
%   function is used to find the operators for more than one enzyme
%   catalyzing the same reaction
% output:
%   selectedOperators: a cummulitave list of all operators to be considered
% function selectedOperators = getEnzymeOperators(enzymeList, rxnID, selectedOperators)
function selectedOperators = getEnzymeOperators(enzymeList, selectedOperators)

global operators
operatorsNum = length(operators);

% if a reaction has more than one enzyme associated to it, operators are
% appended in one list, and below is to keep track of the number of items
% in the operators list so far
if isempty(fieldnames(selectedOperators))
    selectedOperatorsNum = 0;
else
    selectedOperatorsNum = length(selectedOperators);
end

% get the operators of the EC number passed to the function from the
% operators dataset, and operator has to include the enzyme and the
% reaction number to be selected. Some operators don't include all the
% reactions related to them in Kegg
% for enzymeIdx = 1:length(enzymeList)
    for operatorIdx = 1:operatorsNum
%         enzymeMatch = find(strcmp(enzymeList{enzymeIdx}, operators(operatorIdx).Enzyme));
        enzymeMatch = find(strcmp(enzymeList, operators(operatorIdx).Enzyme));
%         rxnMatch = find(operators(operatorIdx).Reaction == rxnID);
        
        if ~isempty(enzymeMatch) %&& ~isempty(rxnMatch)
            selectedOperatorsNum = selectedOperatorsNum+1;
            selectedOperators(selectedOperatorsNum).ID = operators(operatorIdx).ID;
            selectedOperators(selectedOperatorsNum).Reaction = operators(operatorIdx).Reaction;
            selectedOperators(selectedOperatorsNum).Enzyme = operators(operatorIdx).Enzyme;
            selectedOperators(selectedOperatorsNum).Reactant = operators(operatorIdx).Reactant;
            selectedOperators(selectedOperatorsNum).Product = operators(operatorIdx).Product;
            selectedOperators(selectedOperatorsNum).KCF = operators(operatorIdx).KCF;
            
            selectedOperators(selectedOperatorsNum).subKEGGID = operators(operatorIdx).SubstrateID;
            selectedOperators(selectedOperatorsNum).prodKEGGID = operators(operatorIdx).ProductID;
%             selectedOperators(selectedOperatorsNum).CoupledRxn = rxnID;   
        end
    end
% end
end

function [phaseI, phaseImainAtom, phaseImainAtomPos] = sortRs(cypData)
% This function sorts the sturuct based on phase.Reactant.R
for i = 1:length(cypData)
    main{i} = cypData(i).Reactant.R;
end
[sortedR, order] = sort(main);
phaseI = cypData(order);
phaseImainAtom = unique(sortedR);
for i = 1:length(phaseImainAtom)
    for j = 1:length(sortedR)
        if (findstr(phaseImainAtom{i}, sortedR{j}))
            phaseImainAtomPos(i) = j;
            break;
        end
    end
end
end

