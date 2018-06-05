function [isEqual] = compare_two_kcfFiles(kcfFile_1,kcfFile_2)

% function [isEqual] = compare_two_kcfFiles(molFile_1,molFile_2)
% generate KCF files for mol files
% kcfFile_1 = generate_kcfFile(molFile_1);
% kcfFile_2 = generate_kcfFile(molFile_2);
% make sure kcfFiles exist
if ~isempty(kcfFile_1) && ~isempty(kcfFile_2)
    % extract number of atoms and bonds of two mol files
    [atomNum_1,bondNum_1] = extract_atom_bond_number(kcfFile_1);
    [atomNum_2,bondNum_2] = extract_atom_bond_number(kcfFile_2);
    % continue comparing if number of atoms and bonds match
    % and if molecules have more than on atom and have bonds between atoms
    if (atomNum_1 == atomNum_2) && (bondNum_1 == bondNum_2) && (bondNum_1 > 0)
        % extract pairs of atoms and number of bonds between them
        [kcfFile_1_info] = extract_atoms_bonds(kcfFile_1,atomNum_1,bondNum_1);
        [kcfFile_2_info] = extract_atoms_bonds(kcfFile_2,atomNum_2,bondNum_2);
        % compare the mol files information
        [samePairs] = compare_kcfFile_infos(kcfFile_1_info, kcfFile_2_info);
        % if the pairs all match, then molfiles are equal
        if(samePairs)
            isEqual = true;
        else
            isEqual = false;
        end
        % if the number of atoms and bonds do not match, do not compare the details
    else
        isEqual = false;
    end
    % if any of kcf files is empty, there is nothing to compare
else
    isEqual = false;
end
end

function [kcfFile] = generate_kcfFile(molFile)
if isempty(findstr(molFile, '.mol'))
    command = ['curl -F molfile=@', molFile, '.mol -s http://rest.genome.jp/mol2kcf/'];
else
    command = ['curl -F molfile=@', molFile, ' -s http://rest.genome.jp/mol2kcf/'];
end
[status, cmdout] = system(command);
if status ~= 0
    fprintf('Error running curl \n');
    kcfFile = [];
else
    kcfFile = cmdout;
end
end

function [atomNum, bondNum] = extract_atom_bond_number(kcfFile)
% find new lines
newLine_indexes = findstr(kcfFile,sprintf('\n'));
% add index 0 to the array
newLine_indexes = [0 newLine_indexes];
% read line by line
for i = 1 : length(newLine_indexes) - 1
    satrtIndex = newLine_indexes(i) + 1;
    endIndex = newLine_indexes(i + 1) - 1;
    lineText = kcfFile(satrtIndex : endIndex);
    % find number of atoms
    if findstr(lineText,'ATOM')
        atomNum = str2num(lineText(5 : length(lineText)));
    end
    % find number of bonds
    if findstr(lineText,'BOND')
        bondNum = str2num(lineText(5 : length(lineText)));
    end
end
end

function [kcfFile_info] = extract_atoms_bonds(kcfFile,atomNum,bondNum)
% find new lines
newLine_indexes = findstr(kcfFile,sprintf('\n'));
% add index 0 to the array
newLine_indexes = [0 newLine_indexes];
% read line by line
for i = 1 : length(newLine_indexes) - 1
    satrtIndex = newLine_indexes(i) + 1;
    endIndex = newLine_indexes(i + 1) - 1;
    lineText = kcfFile(satrtIndex : endIndex);
    % find where to start reading atom types
    if findstr(lineText,'ATOM')
        atomLines_count = 1;
        readAtomTypes_startIndex = i + 1;
    end
    % find where to start reading bond numbers
    if findstr(lineText,'BOND')
        BondLines_count = 1;
        readBonds_startIndex = i + 1;
    end
end
% read atom types
for i = readAtomTypes_startIndex : readAtomTypes_startIndex + (atomNum - 1)
    % read a line
    satrtIndex = newLine_indexes(i) + 1;
    endIndex = newLine_indexes(i + 1) - 1;
    lineText = kcfFile(satrtIndex : endIndex);
    % split the content of the line
    lineContent = strsplit(lineText,' ');
    % remove empty cells
    lineContent_new = lineContent(~cellfun(@isempty, lineContent));
    % the forth cell is always the type of the atom
    atomTypes(atomLines_count) = lineContent_new(2);
    atomLines_count = atomLines_count + 1;
end
% generate the pairs
for i = readBonds_startIndex : readBonds_startIndex + (bondNum - 1)
    % read a line
    satrtIndex = newLine_indexes(i) + 1;
    endIndex = newLine_indexes(i + 1) - 1;
    lineText = kcfFile(satrtIndex : endIndex);
    % split the content of the line
    lineContent = strsplit(lineText,' ');
    % remove empty cells
    lineContent_new = lineContent(~cellfun(@isempty, lineContent));
    % extract the pair information
    atom1 = atomTypes(str2num(lineContent_new{2}));
    atom2 = atomTypes(str2num(lineContent_new{3}));
    bond = str2num(lineContent_new{4});
    % store the information
    kcfFile_info(BondLines_count).atom1 = atom1;
    kcfFile_info(BondLines_count).atom2 = atom2;
    kcfFile_info(BondLines_count).bond = bond;
    BondLines_count = BondLines_count + 1;
end
end

function [samePairs] = compare_kcfFile_infos(kcfFile_1_info, kcfFile_2_info)
% initialization
foundMatch = zeros(1,length(kcfFile_1_info));
alreadyMatched = [];
% check the apirs
for i = 1 : length(kcfFile_1_info)
    found = false;
    atom1_f1 = kcfFile_1_info(i).atom1;
    atom2_f1 = kcfFile_1_info(i).atom2;
    bond_f1 = kcfFile_1_info(i).bond;
    for j = 1 : length(kcfFile_2_info)
        % check if the row has not been matched before
        if isempty(find(alreadyMatched == j, 1)) && (found == false)
            atom1_f2 = kcfFile_2_info(j).atom1;
            atom2_f2 = kcfFile_2_info(j).atom2;
            bond_f2 = kcfFile_2_info(j).bond;
            % check atom types
            if (strcmp(atom1_f1,atom1_f2) && (strcmp(atom2_f1,atom2_f2)))||(strcmp(atom1_f1,atom2_f2) && (strcmp(atom2_f1,atom1_f2)))
                % check bond numbers
                if (bond_f1 == bond_f2)
                    found = true;
                    foundMatch(i) = 1;
                    alreadyMatched = [alreadyMatched j];
                end
            end
        end
    end
end
% check if all the pairs matched to something
if (prod(foundMatch) == 1)
    samePairs = true;
else
    samePairs = false;
end
end