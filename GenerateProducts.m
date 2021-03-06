% Author: Mona YousoufShahi

function [inputList, inputStructure, inputListNumbering, products, operatorIdx] = GenerateProducts(inputMolecule, operators)
% inputMolecule: The kegg ID for the new compound
% operators: filename, contains operator's information
% inputList: input atom list
% products: outputs after applying operators

% construct a list from the input
% first column is main atomtype and others are neighbors
% inputMolecule = 'C00463';    %kegg id for indole
[inputList, inputStructure, inputListNumbering] = ConstructInputListfromKEGG(inputMolecule);

% for each main and its neighbors, output the possible products
[products, operatorIdx] = ConstructOutputList(inputList, inputStructure, operators);
end

function [inputList, inputStructure, inputListNumbering] = ConstructInputListfromKEGG(inputMolecule)
inputMolecule

if strcmp(lower(inputMolecule(1)), 'c') && (length(inputMolecule) == 6)
    inputStructure = GetKCF(inputMolecule);
else
    inputStructure = ReadMolFile(inputMolecule);
end
inputStructure.Atoms(:,2) = Mutate(inputStructure.Atoms(:,2));
[inputList, inputListNumbering]= FindAllRMs(inputStructure);
end


function s = Mutate(s)
    if isstruct(s)
        if length(s) == 1
            fields = fieldnames(s);
            for i = 1:length(fields)
                field = fields{i};
                s.(field) = Mutate(s.(field));
            end
        else
            for i = 1:length(s)
                s(i) = Mutate(s(i));
            end
        end
    elseif iscell(s)
        for i = 1:length(s)
            s{i} = Mutate(s{i});
        end
    else
        keggatomTypeOld = {'C2x', 'C2y', 'N4x', 'N2x', 'N1y'};
        keggatomTypeNew = {'C8x', 'C8y', 'N1x', 'N5x', 'N4y'};
        idx = find(strcmp(s, keggatomTypeOld) == 1);
        if ~isempty(idx)
            s = keggatomTypeNew{idx};
        end
    end
end

function molecule = GetKCF(inputMolecule)
% retrieving the mol file from kegg databse
url = ['http://rest.kegg.jp/get/', inputMolecule, '/kcf'];
str = urlread(url);
molecule = ReadKCF(str);
end

function [products, operatorIdx] = ConstructOutputList(inputList, inputStructure, operators)
load(operators)

% applying operators
[products, operatorIdx] = ApplyOperators(inputList, selectedOperators, operatorsMainAtom, operatorsMainAtomPos);
% we can use inputStructure to construct an input list from products
end

function [op, operatorIdx] = ApplyOperators(inputList, operators, operatorsMainAtom, operatorsMainAtomPos)
nAtoms = length(inputList);
counterNewProd = 1;

operatorIdx = [];
for i = 1:nAtoms
    % comparing the main atom of input with the main atoms of phaseI
    indexAT = find(ismember(operatorsMainAtom, inputList(i).R)==1);
    if isempty(indexAT)
        continue;
    end
    
    if length(indexAT) == 1 && indexAT == length(operatorsMainAtom)
        endPattern = operatorsMainAtomPos(indexAT);
    else
        endPattern = operatorsMainAtomPos(indexAT+1)-1;
    end
    
    searchEntries = operators(operatorsMainAtomPos(indexAT):endPattern);
    searchEntriesIdx = operatorsMainAtomPos(indexAT):endPattern;
    for j = 1:length(searchEntries)
        matchedPattern = true;
        for k = 1:4
            comp1 = eval(['searchEntries(j).Reactant.M', num2str(k)]);
            comp2 = eval(['inputList(i).M', num2str(k)]);
            if isempty(comp1) && isempty(comp2) && matchedPattern
                matchedPattern = true;
            elseif strcmp(comp1, comp2) && matchedPattern
                matchedPattern = true;
            else
                matchedPattern = false;
                break
            end
            
%  Comment out this section, so the tool would only check for the neighbors
%  of the center metabolite, and not check for neighbors of beighbors

%             comp1 = eval(['searchEntries(j).Reactant.M', num2str(k), 'N']);
%             comp2 = eval(['inputList(i).M', num2str(k), 'N']);
%             if isempty(comp1) && isempty(comp2) && matchedPattern
%                 matchedPattern = true;
%             elseif isequal(comp1, comp2) && matchedPattern
%                 matchedPattern = true;
%             else
%                 matchedPattern = false;
%                 break
%             end
        end
        if matchedPattern
            searchEntries(j).ObtainedFromInput = i;
            op(counterNewProd) = searchEntries(j);
            operatorIdx(counterNewProd) = searchEntriesIdx(j);
            counterNewProd = counterNewProd + 1;
        end
    end
end
if counterNewProd == 1
    op = [];
end
end

function molecule = ReadMolFile(filenames)
% reads mol files and extracts both atoms and bonds
if isempty(findstr(filenames, '.mol'))
    command = ['curl -F molfile=@', filenames, '.mol -s http://rest.genome.jp/mol2kcf/'];
else
    command = ['curl -F molfile=@', filenames, ' -s http://rest.genome.jp/mol2kcf/'];
end
[status, cmdout] = system(command);
if status ~= 0
    fprintf('Error running curl \n');
else
    molecule = ReadKCF(cmdout);
end
end

function molecule = ReadKCF(cmdout)
    sIndexAtom = findstr(cmdout, 'ATOM');
    eIndexAtom = findstr(cmdout, 'BOND');
    atomText = cmdout(sIndexAtom(1)+length('ATOM'):eIndexAtom(1)-1);
    atomLen = sscanf(atomText, '%d', 1);
    molecule.Atoms = ReadAtom(atomText, atomLen);
    bondText = cmdout(eIndexAtom(1)+length('BOND'):end);
    bondLen = sscanf(bondText, '%d', 1);
    if bondLen ~= 0
        molecule.Bonds = ReadBond(bondText, bondLen);
    else
        molecule.Bonds = [];
    end
end

function atom = ReadAtom(text, len)
indexNewLineAlign = find(text == 10);
for atomCount = 1:len
    oneLine = text(indexNewLineAlign(atomCount)+1:indexNewLineAlign(atomCount+1));
    components = sscanf(oneLine, '%d %s %s %f %f');
    indexNotSpace = find(oneLine ~= 32);
    repeatedPos = find((diff(indexNotSpace)==1)==0);
    lenChar1 = repeatedPos(2)-repeatedPos(1);
    lenChar2 = repeatedPos(3)-repeatedPos(2);
    atom{atomCount, 1} = components(1);
    atom{atomCount, 2} = char(components(2:2+lenChar1-1))';
    atom{atomCount, 3} = char(components(2+lenChar1:2+lenChar1+lenChar2-1))';
    atom{atomCount, 4} = components(end-1);
    atom{atomCount, 5} = components(end);
end
end

function bond = ReadBond(text, len)
indexNewLineAlign = find(text == 10);
for atomCount = 1:len
    oneLine = text(indexNewLineAlign(atomCount)+1:indexNewLineAlign(atomCount+1));
    components = sscanf(oneLine, '%d %d %d %d');
    bond(atomCount, 1:4) = components;
end
end

function [inputList inputListNumbering] = FindAllRMs(inputStructure)
% For each atom, finds the neighbors and put them in a row of inputList

% Bond: 1st&2nd cols-atom#, 3rd number of bonds
for i = 1:size(inputStructure.Atoms, 1)
    % first col of inputList is main atom
    inputList(i).R = inputStructure.Atoms{i, 2};
    [nghrs tempVar] = FindNeighbors(inputStructure.Bonds, i);
    neighbors = inputStructure.Atoms(nghrs, 2)';
    for j = 1:length(neighbors)
        eval(['inputList(i).M', num2str(j), ' = neighbors{j};']);
        [nghrsN tempVar] = FindNeighbors(inputStructure.Bonds, nghrs(j));
        neighborsN = inputStructure.Atoms(nghrsN, 2)';
        neighborsN = sort(neighborsN);
        eval(['inputList(i).M', num2str(j), 'N = neighborsN;']);
    end
    for j = length(neighbors)+1:4
        eval(['inputList(i).M', num2str(j), ' = [];']);
        eval(['inputList(i).M', num2str(j), 'N = [];']);
    end    
end
    [inputList, inputListNumbering] = SortFields(inputList, inputStructure);
end

function [output, inputListNumbering] = SortFields(input, inputStructure)
inputListNumbering = zeros(size(inputStructure.Atoms, 1), 5);
output = input;
for i = 1:length(input)
    inputListNumbering(i, 1) = i;
%     clear tempData
    tempData = {};
    
    for j = 1:4
        if (eval(['isempty(input(i).M', num2str(j), ')'])) || (eval(['isempty(input(i).M', num2str(j), 'N)']))
            break
        end
        command = ['tempData{j} = strcat(input(i).M', num2str(j), ', input(i).M', num2str(j), 'N{1:end});'];
        eval(command);
    end
    [neighboringAtom neighborsPos] = FindNeighbors(inputStructure.Bonds, i);
    for j = 1:length(neighborsPos)
        tempData{j}(end+1) = num2str(inputStructure.Bonds(neighborsPos(j), 4));
    end
    [sortedMs order] = sort(tempData);
    inputListNumbering(i, 2:length(order)+1) = neighboringAtom(order);
    for j = 1:length(order)
        eval(['output(i).M', num2str(j),' = input(i).M', num2str(order(j)), ';']);
        eval(['output(i).M', num2str(j),'N = input(i).M', num2str(order(j)), 'N;']);
    end
end
end

function [neighboringAtom neighborsPos] = FindNeighbors(bonds, atom)
% Finds the neighbors of an atom
col2 = find(bonds(:, 2) == atom);
col3 = find(bonds(:, 3) == atom);
neighborsPos = unique([col2; col3]);
k = 1; neighboringAtom = [];
for x = neighborsPos.'
    if bonds(x, 3) == atom
        neighboringAtom(k) = bonds(x,2);
    else
        neighboringAtom(k) = bonds(x,3);
    end
    k = k + 1;
end
end
