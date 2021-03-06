% This function retrieves a unique list of EC numbers saved in the model
function enzymeList = getEnzymeList(mainModel)

genesNum = length(mainModel.genes);
count = 1;

% if a gene is associated with more than one EC number, split those EC
% numbers in separate entities
for i = 1:genesNum
    if ~isempty(mainModel.ECNums{i})
         ecStr = strsplit(char(mainModel.ECNums{i}),';');
         ecSize = length(ecStr);
         if  ecSize > 1
             for j = 1:ecSize
                 ecList{count,1} = ecStr{j}; 
                 count = count + 1;
             end
             
         else
             ecList{count,1} = ecStr{1}; 
             count = count + 1;
         end
    end
end

% Remove duplicated EC numbers 
enzymeList = unique(ecList);
end