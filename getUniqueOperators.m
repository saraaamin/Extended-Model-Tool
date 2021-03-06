% This function takes the set of generated operators and remove duplicate
% enteries of the same operator
function uniqueOperators = getUniqueOperators(operators)

uniqueOperators(1) = operators(1);
index = 1;

for i = 2:length(operators)
    duplicate = false;
        similarIDx = isSimilar(uniqueOperators, operators(i));
        if isempty(similarIDx)
            index = index +1;
            uniqueOperators(index) = operators(i);
        end    
end
    
end


% This function checks of the details of two operators are the same so they
% would be considered duplicates
function similarIDx = isSimilar(all_operators, queryOperator)
similarIDx = [];

for i = 1:length(all_operators)
    % if the length of RPairs IDs are equal then there is a chance that the
    % current operator and the query operator are the same
    if length(all_operators(i).ID) == length(queryOperator.ID)
        
        % if both IDs are equal, check if other fields are equal
        if isempty(find((all_operators(i).ID == queryOperator.ID) == 0))
            if isempty(find((all_operators(i).subKEGGID == queryOperator.subKEGGID) == 0))
                if isempty(find((all_operators(i).prodKEGGID == queryOperator.prodKEGGID) == 0))
                    similarIDx = i;
                    return;
                end
            end
        end
    end            
end
end
        