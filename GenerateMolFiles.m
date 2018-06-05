%  A function to generate a Mol file for a specific product
function [prodEnzymesList, operatorIdx] = GenerateMolFiles(inputStructure, inputListNumbering, products, opsIdx)

% idx = strfind(targetMolecule, '.mol');
% if isempty(idx)
%     mkdir(targetMolecule);
% else
%     targetMolecule = targetMolecule(1:idx-1);
%     mkdir(targetMolecule);
% end
targetMolecule = 'productsMol';
mkdir(targetMolecule);
delete([targetMolecule, '\*.mol']);
delete([targetMolecule, '\*.smiles']);

reactionListRemoved = [74 4510 3472 7066 7079];

prodEnzymesList = [];
index = 1;
for j = 1:length(products)
    if isempty(products(j).KCF.M) || ~isempty(intersect(products(j).Reaction, reactionListRemoved))
        continue;
    end
    % adding the added structure to inputStructure
    modifiedInputStructure = BuildNewStructure(products(j), inputStructure, inputListNumbering(products(j).ObtainedFromInput, :));
    % constructing the mol file
    ConstructMolFile(modifiedInputStructure, [targetMolecule, '\product_', num2str(index), '.mol']); % windows
%     ConstructMolFile(modifiedInputStructure, [targetMolecule, '/product_', num2str(j), '.mol']);  % unix
    prodEnzymesList{index,1} = products(j).Enzyme;
    operatorIdx(index,1) = opsIdx(j);
    index = index + 1;
end
end

function modifiedInputStructure = BuildNewStructure(pattern, inputStructure, inputNumbering)
% pattern: pattern applied to R atom
% inputStructure: atom and bond information for the test molecule
% inputNumbering: array with first element is R, others are Ms

err = 1;
modifiedInputStructure = inputStructure;

% changing the atom type for R
modifiedInputStructure = ChangeR(modifiedInputStructure, pattern, inputNumbering);
if ~isempty(find(modifiedInputStructure.Bonds==0))
    err = 0;
end

% updaing the bonds between R-Ms and Ms-neighbors
modifiedInputStructure = UpdateBonds(modifiedInputStructure, pattern, inputNumbering);
if ~isempty(find(modifiedInputStructure.Bonds==0))
    err = 0;
end

% changing Ms
modifiedInputStructure = ChangeMs(modifiedInputStructure, pattern, inputNumbering);
if ~isempty(find(modifiedInputStructure.Bonds==0))
    err = 0;
end

% changing D
modifiedInputStructure = ChangeD(modifiedInputStructure, pattern, inputNumbering, inputStructure);
if ~isempty(find(modifiedInputStructure.Bonds==0))
    err = 0;
end

if err == 0 || isempty(modifiedInputStructure.Bonds)
    return
end

% removing atoms not connected to the structure
modifiedInputStructure = RemovePartsOfStructure(modifiedInputStructure, pattern.ObtainedFromInput);
end

function oStruct = ChangeR(oStruct, pattern, RMList)
% change the atom type for R
oStruct.Atoms{RMList(1), 2} = pattern.Product.R;
end

function oStruct = UpdateBonds(oStruct, pattern, RMList)
patternRMbonds = FindRMnumberofBonds(pattern.KCF);
% change the number of bonds between R and Ms if the number of neighbors
% are the same
if sum(RMList~=0)-1 == sum(patternRMbonds~=0)
    sumB = 0;
    for i = 2:5
        if (RMList(i)~=0)
            indexRM(i-1) = FindBondBetween2atoms(oStruct.Bonds, RMList(1), RMList(i));
            sumB = sumB + oStruct.Bonds(indexRM(i-1), 4);
        end
    end
    if (sum(patternRMbonds) ~= sumB)
        for i = 1:length(indexRM)
            oStruct.Bonds(indexRM(i), 4) = patternRMbonds(i);
        end
        % find the atoms in the loop (1)R, (2)M1, M2, (3)N1, N2, (4)C
        %    R
        % M1   M2
        % N1   N2
        %    C
        [cycle atoms bFound] = FindCycle(oStruct.Bonds, RMList(1));
        
        if bFound && length(atoms) == 7 && ~isempty(cycle == RMList(1))% checks if the cycle generates a hexagonal and the cycle includes R
            % if R-M1 and R-M2 are 1, M1-N1 and M2-N2 should be 2
            % N1-C and N2-C should be 1
            if oStruct.Bonds(cycle(1), 4)== 1 && oStruct.Bonds(cycle(6), 4) == 1
                oStruct.Bonds(cycle(2), 4) = 2;
                oStruct.Bonds(cycle(5), 4) = 2;
                oStruct.Bonds(cycle(3), 4) = 1;
                oStruct.Bonds(cycle(4), 4) = 1;
            end
            % if R-M1 is 1 and R-M2 is 2, then M1-N1 is 2, M2-N2 is 1, N2-C is 2, N1-C
            % is 1
            if oStruct.Bonds(cycle(1), 4)== 2 && oStruct.Bonds(cycle(6), 4) == 1
                oStruct.Bonds(cycle(3), 4) = 2;
                oStruct.Bonds(cycle(5), 4) = 2;
                oStruct.Bonds(cycle(2), 4) = 1;
                oStruct.Bonds(cycle(4), 4) = 1;
            end
            
            % if R-M1 is 2 and R-M2 is 1, then M1-N1 is 1, M2-N2 is 2, N2-C is 1, N1-C
            % is 2
            if oStruct.Bonds(cycle(1), 4)== 1 && oStruct.Bonds(cycle(6), 4) == 2
                oStruct.Bonds(cycle(2), 4) = 2;
                oStruct.Bonds(cycle(4), 4) = 2;
                oStruct.Bonds(cycle(3), 4) = 1;
                oStruct.Bonds(cycle(5), 4) = 1;
            end
        end
    end
% by Neda: update the bonds between R and Ms even if the number of neighbors are not the same    
else
    indexRM = [];
    for i = 2:5
        if (RMList(i)~=0) && (patternRMbonds(i-1)~=0)
            indexRM(i-1) = FindBondBetween2atoms(oStruct.Bonds, RMList(1), RMList(i));
            
        end
    end
    if ~isempty(indexRM)
        for i = 1:length(indexRM)
            if (indexRM(i) ~= 0)
                oStruct.Bonds(indexRM(i), 4) = patternRMbonds(i);
            end
        end
    end
end
% finish by Neda
end

function oStruct = ChangeMs(oStruct, pattern, RMList)

% removing the bond between a removed neighbor and R
for z = 2:length(RMList)
    if RMList(z) ~= 0
        if eval(['pattern.Product.M', num2str(z-1), '== ''*'''])  % if the neighbor changes to *, its bond to R should be removed
            indexRM1 = FindBondBetween2atoms(oStruct.Bonds, RMList(1), RMList(z));
            oStruct.Bonds(indexRM1, :) = [];
            nBonds = size(oStruct.Bonds, 1);
            oStruct = UpdateBondsM(oStruct, nBonds, indexRM1);
            % removing the atom if it is not connected to other atoms
            Mneighbor1 = FindNeighbors(oStruct.Bonds, RMList(z));
            if isempty(Mneighbor1)
                oStruct = RemoveAtom(oStruct, RMList(z));
                % update the RMList
                for k = z+1:length(RMList)
                    if RMList(k) ~= 0 && RMList(k) > RMList(z)
                        RMList(k) = RMList(k)-1;
                    end
                end
            end
            % some cases a group of atoms need to be removed
            % any atom not connected to R should be removed
            
        end
    end
end
end

function oStruct = UpdateBondsM(oStruct, nBonds, indexBond)
for m = indexBond:nBonds
    oStruct.Bonds(m, 1) = m;
end
end

function oStruct = RemoveAtom(oStruct, atomTobeRemoved)
% remove atoms
oStruct.Atoms(atomTobeRemoved, :) = [];
% update atom numbers
for m = atomTobeRemoved:size(oStruct.Atoms, 1)
    oStruct.Atoms{m, 1} = m;
end
% remove bonds
for m = 1:length(atomTobeRemoved)
    index = FindNeighborsP(oStruct.Bonds, atomTobeRemoved(m));
    oStruct.Bonds(index, :) = [];
end
% update atoms in bonds
for i = 1:size(oStruct.Bonds, 1)
    oStruct.Bonds(i, 1) = i;
    oStruct.Bonds(i, 2) = oStruct.Bonds(i, 2) - sum((atomTobeRemoved<oStruct.Bonds(i, 2)));
    oStruct.Bonds(i, 3) = oStruct.Bonds(i, 3) - sum((atomTobeRemoved<oStruct.Bonds(i, 3)));
end
end

function oStruct = ChangeD(oStruct, pattern, RMList, iStruct)
oStruct = AddNewFunctionalGroups(oStruct, pattern, RMList, iStruct);
% % % in some cases, Ms can have new functional groups
% allMList = FindNeighbors(oStruct.Bonds, RMList(1));
% MList = RMList(2:end);
% MListdiff = setdiff(MList(MList(2:end)~=0), allMList);
% for i = 1:length(MListdiff)
%     for j = 1:4
%         if strcmp(eval(['pattern.Product.M',num2str(j)]), oStruct.Atoms(MList(i), 2))
%             neighbors = FindNeighbors(oStruct.Bonds, MList(i));
%             pneighbors = eval(['length(pattern.Product.M', num2str(j), 'N);']);
%             if (length(neighbors) ~= pneighbors)
%                 test = 1;
%             end
%         end
%     end
% end

% checking whether the number of distant neighbor is changed
for i = 1:4
    if (RMList(i+1) ~= 0)
        if eval(['length(pattern.Product.M', num2str(i), 'N) > length(pattern.Reactant.M', num2str(i), 'N)'])
            MNList = zeros(1, 5);
            MNList(1) = RMList(i+1);
            Mneighbors = FindNeighbors(oStruct.Bonds, MNList(1));
            MNList(2:length(Mneighbors)+1) = Mneighbors;
            if eval(['length(pattern.Product.M', num2str(i), 'N) ~= length(Mneighbors)'])
                % current neighbors for M_i
                currNeighbors = eval(['pattern.Product.M', num2str(i), 'N']);
                % find new neighbors for M_i
                atomD = setdiff(currNeighbors, oStruct.Atoms(Mneighbors, 2)');
                atomD(find(strcmp(atomD, '*')==1)) = [];
                % new D's index
                newDindex1 = find(ismember(pattern.KCF.atom2(:, 2)', atomD) == 1);
                newDindex2 = FindNeighbors(pattern.KCF.bond2, pattern.KCF.M(i, 2));
                indexMN = intersect(newDindex1, newDindex2);
                if isempty(indexMN)
                    continue
                end
                if ~isempty(pattern.KCF.D)
                    for z = 1:size(pattern.KCF.D, 1)
                        indexMN(find(pattern.KCF.D(z, 2)==indexMN))= [];
                    end
                end
                % remove atoms that can be aligned with rectant atoms (finding D atoms)
                if length(indexMN)>1
                    for a = length(indexMN):-1:1
                        idx = find(indexMN(a) == pattern.KCF.Align(:, 2));
                        if ~isempty(idx) & pattern.KCF.Align(idx, 1) ~= -1
                            indexMN(a) = [];
                        end
                    end
                end
                patternNew.Product.D1.D = atomD;
                % error happens here if ...
                if length(indexMN)>1 || length(atomD)>1 || isempty(indexMN)
                    continue;
                end
                patternNew.Product.D1.atomD = indexMN;
                % find atom R
                atomR = FindNeighbors(pattern.KCF.bond2, indexMN);
                % need to handle this case (e.g. look at R00074), how to
                % handle this?
                if length(atomR)>1
                    continue
                end
                % remove atoms that cannot be aligned with rectant atoms (meaning atoms are part of D)
                for a = length(atomR):-1:1
                    if isempty(find(atomR(a) == pattern.KCF.Align(:, 2)))
                        atomR(a) = [];
                    end
                end
                patternNew.Product.D1.atomR = atomR;
                patternNew.Product.D1.atom = pattern.KCF.atom2(indexMN, :);
                patternNew.Product.D1.bond = pattern.KCF.bond2(FindBondBetween2atoms(pattern.KCF.bond2, patternNew.Product.D1.atomD, patternNew.Product.D1.atomR), :);
                patternNew.Product.D1.bond(1) = 1;
                oStruct = AddNewFunctionalGroups(oStruct, patternNew, MNList, iStruct);
            end
        end
    end
end

end

function oStruct = RemovePartsOfStructure(oStruct, nodeR)
% sparse adjacency matrix
A = oStruct.Bonds;
adjMatrix = sparse([A(:,2);A(:,3)], [A(:,3);A(:,2)], 1, size(oStruct.Atoms, 1), size(oStruct.Atoms, 1));
connectedNodes = graphtraverse(adjMatrix, nodeR);
unConnectedNodes = setdiff(1:size(oStruct.Atoms,1), connectedNodes);
oStruct = RemoveAtom(oStruct, unConnectedNodes);
end

function oStruct = AddNewFunctionalGroups(oStruct, pattern, RMList, iStruct)
nBonds = size(oStruct.Bonds, 1);
nAtoms = size(oStruct.Atoms, 1);
if ~isempty(pattern.Product.D1.D) && ~isempty(pattern.Product.D1.atom)
    A = 0; B = 0; C = 0; D = 0;
    O = [iStruct.Atoms{RMList(1), 4} iStruct.Atoms{RMList(1), 5}];
    % finds O2x position
    if strcmp(pattern.Product.D1.D, 'O2x')
        if sum(strcmp(pattern.Product.M1N, 'O2x'))>0 || RMList(3)==0
            A = [iStruct.Atoms{RMList(2), 4} iStruct.Atoms{RMList(2), 5}];
        else
            A = [iStruct.Atoms{RMList(3), 4} iStruct.Atoms{RMList(3), 5}];
        end
    else
        A = [iStruct.Atoms{RMList(2), 4} iStruct.Atoms{RMList(2), 5}];
        if RMList(3) ~= 0
            B = [iStruct.Atoms{RMList(3), 4} iStruct.Atoms{RMList(3), 5}];
        end
        if RMList(4) ~= 0
            C = [iStruct.Atoms{RMList(4), 4} iStruct.Atoms{RMList(4), 5}];
        end
        if RMList(5) ~= 0
            D = [iStruct.Atoms{RMList(5), 4} iStruct.Atoms{RMList(5), 5}];
        end
    end
    lenX = 0.8;
    % Finding the position of D
    X = getThirdVector(O, lenX, A, B, C, D);
    % making larger distance between R and D
    %    X = X + 0.7*(X-O);
    %    X = X + 0.556*(X-O);
    for i = 1:size(pattern.Product.D1.atom, 1)
        atomList(i) = pattern.Product.D1.atom{i, 1};
    end
    atomD = find(atomList == pattern.Product.D1.atomD);
    
    % adding new atoms (D)
    atom = nAtoms+1;
    atomMapping = [pattern.Product.D1.atomR RMList(1)];
    for i = 1:size(pattern.Product.D1.atom, 1)
        atomMapping = [atomMapping; pattern.Product.D1.atom{i, 1} atom];
        oStruct.Atoms{atom, 1} = atom;
        oStruct.Atoms(atom, 2:3) = pattern.Product.D1.atom(i, 2:3);
        pos1 = pattern.Product.D1.atom{i, 4} - pattern.Product.D1.atom{atomD, 4};
        pos2 = pattern.Product.D1.atom{i, 5} - pattern.Product.D1.atom{atomD, 5};
        oStruct.Atoms{atom, 4} = X(1) + pos1;
        oStruct.Atoms{atom, 5} = X(2) + pos2;
        atom = atom + 1;
    end
    
    % rotate the functional group based on position of R, D and D's neighbors
    atomDneighbor = FindNeighbors(pattern.Product.D1.bond, pattern.Product.D1.atomD);
    if ~isempty(atomDneighbor)
        atomDneighbor(find(atomDneighbor == pattern.Product.D1.atomR)) = [];
    end
    if ~isempty(atomDneighbor) && ~isempty(find(atomMapping(:, 1) == atomDneighbor(1)))
        atomDneighborinoStruct = atomMapping(find(atomMapping(:, 1) == atomDneighbor(1)), 2);
        p0 = X; % D
        p1 = O; % R
        p2(1) = oStruct.Atoms{atomDneighborinoStruct, 4};   % D's neighbor in the added functional group
        p2(2) = oStruct.Atoms{atomDneighborinoStruct, 5};
        % find if functional group is at right/left of D
        indexAtomMapping = atomMapping(2:end, 2);
        indexAtomMapping(find(indexAtomMapping==RMList(1), 1)) = [];
        minX = min([oStruct.Atoms{atomMapping(2:end, 2), 4}]);
        maxX = max([oStruct.Atoms{atomMapping(2:end, 2), 4}]);
        bTheta = true;
        if abs(minX-p0(1))<abs(maxX-p0(1))   % D[....]
            if p0(1)<p1(1)  % D--R
                bTheta = false;
            else            % R--D
                bTheta = true;
            end
        elseif abs(minX-p0(1))>abs(maxX-p0(1)) % [....]D
            if p0(1)>p1(1)  % R--D
                bTheta = false;
            else            % D--R
                bTheta = true;
            end
        end
        if ~bTheta
            count = 1; atomsChanged = [];
            for i = 1:size(atomMapping, 1)
                if atomMapping(i, 1) == pattern.Product.D1.atomD || atomMapping(i, 1) == pattern.Product.D1.atomR
                    continue
                end
                xpos(count) = oStruct.Atoms{atomMapping(i, 2), 4};
                ypos(count) = oStruct.Atoms{atomMapping(i, 2), 5};
                atomsChanged = [atomsChanged atomMapping(i, 2)];
                count = count + 1;
            end
            bRotation = false;
            x_rotated = XMirrorGroupofPoints(xpos, p0(1));
            y_rotated = ypos;
            bRotation = true;
            if bRotation
                for i = 1:length(atomsChanged)
                    oStruct.Atoms{atomsChanged(i), 4} = x_rotated(i);
                    oStruct.Atoms{atomsChanged(i), 5} = y_rotated(i);
                end
            end
        end
        
        %         if bTheta
        %             theta = 0;
        %         else
        %             theta = 180 - atan2(abs(det([p2-p0;p1-p0])),dot(p2-p0,p1-p0))*180/pi;
        %         end
        %         if abs(theta) > 10
        %             count = 1; atomsChanged = [];
        %             for i = 1:size(atomMapping, 1)
        %                 if atomMapping(i, 1) == pattern.Product.D1.atomD || atomMapping(i, 1) == pattern.Product.D1.atomR
        %                     continue
        %                 end
        %                 xpos(count) = oStruct.Atoms{atomMapping(i, 2), 4};
        %                 ypos(count) = oStruct.Atoms{atomMapping(i, 2), 5};
        %                 atomsChanged = [atomsChanged atomMapping(i, 2)];
        %                 count = count + 1;
        %             end
        %             bRotation = false;
        %             [x_rotated y_rotated] = RotateGroupofPoints(xpos, ypos, p0(1), p0(2), theta);
        %             bRotation = true;
        %             if bRotation
        %                 for i = 1:length(atomsChanged)
        %                     oStruct.Atoms{atomsChanged(i), 4} = x_rotated(i);
        %                     oStruct.Atoms{atomsChanged(i), 5} = y_rotated(i);
        %                 end
        %             end
        %         end
    end
    
    % adding the bond between R and D
    bondNum = nBonds+1;
    oStruct.Bonds(bondNum, 1) = bondNum;
    oStruct.Bonds(bondNum, 2) = RMList(1);
    indexD = find(atomMapping(:, 1) == pattern.Product.D1.atomD);
    oStruct.Bonds(bondNum, 3) = atomMapping(indexD, 2);
    RDRowinBond = FindBondBetween2atoms(pattern.Product.D1.bond, pattern.Product.D1.atomR, pattern.Product.D1.atomD);
    oStruct.Bonds(bondNum, 4) = pattern.Product.D1.bond(RDRowinBond, 4);
    bondNum = bondNum + 1;
    
    % adding other bonds in the functional group D
    for i = 1:size(pattern.Product.D1.bond, 1)
        bRemoveBond = false;
        if i ~= RDRowinBond
            oStruct.Bonds(bondNum, 1) = bondNum;
            for j = 2:3
                indexAtom = find(atomMapping(:, 1) == pattern.Product.D1.bond(i, j));
                if ~isempty(indexAtom)
                    oStruct.Bonds(bondNum, j) = atomMapping(indexAtom, 2);
                elseif strcmp(pattern.Product.D1.D, 'O2x')       % case of o2x where the atom is connect to two other atoms
                    neighborsM1 = FindNeighbors(iStruct.Bonds, RMList(2));
                    neighborsM2 = FindNeighbors(iStruct.Bonds, RMList(3));
                    connectToM = 0;
                    if length(pattern.Reactant.M1N) < length(pattern.Product.M1N)
                        % M1 is connected to O2x
                        if isequal(sort(iStruct.Atoms(neighborsM1, 2)'), sort(pattern.Reactant.M1N))
                            connectToM = RMList(2);
                        elseif isequal(sort(iStruct.Atoms(neighborsM2, 2)'), sort(pattern.Reactant.M1N))
                            connectToM = RMList(3);
                        end
                    elseif length(pattern.Reactant.M2N) < length(pattern.Product.M2N)
                        % M2 is connected to O2x
                        if isequal(sort(iStruct.Atoms(neighborsM1, 2)'), sort(pattern.Reactant.M2N))
                            connectToM = RMList(2);
                        elseif isequal(sort(iStruct.Atoms(neighborsM2, 2)'), sort(pattern.Reactant.M2N))
                            connectToM = RMList(3);
                        end
                    else
                        bRemoveBond = true;
                    end
                    oStruct.Bonds(bondNum, j) = connectToM;
                else
                    bRemoveBond = true;
                end
            end
            oStruct.Bonds(bondNum, 4) = pattern.Product.D1.bond(i, 4);
            if ~bRemoveBond
                bondNum = bondNum +1;
            end
        end
    end
end
end

function xpos = XMirrorGroupofPoints(xpos, x_center)
xpos = xpos-2*(xpos-x_center);
end

function X = getThirdVector(O, lenX, A, B, C, D)
if A~=0 & B==0 & C==0 & D==0
    v = [O(1) A(1); O(2) A(2)];
    x_center = O(1);
    y_center = O(2);
    center = repmat([x_center; y_center], 1, 2);
    theta = 2*pi/3;
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    vo = R*(v - center) + center;
    X(1) = vo(1,2);
    X(2) = vo(2,2);
else
    nA = A - O;
    nA = nA / norm(nA);
    nB = B - O;
    nB = nB / norm(nB);
    nC = C - O;
    nC = nC / norm(nC);
    nD = D - O;
    nD = nD / norm(nD);
    if A == 0, nA=0; end
    if B == 0, nB=0; end
    if C == 0, nC=0; end
    if D == 0, nD=0; end
    nX = -(nA + nB + nC + nD);
    nX = nX / norm(nX);
    X = O + nX * lenX;
end
end

function patternRMbonds = FindRMnumberofBonds(kcf)
% finding the number of bonds between R and Ms in atom2

R = kcf.R(1, 2);
M = zeros(4, 1);
M(1:size(kcf.M, 1)) = kcf.M(:, 2);

for i = 1:4
    if M(i) > 0
        % need to check if the atoms are the same (they are not sorted so we might have it in another place)
        %         if isequal(pattern.KCF.atom1{M(i), 1}, oStruct.Atoms{RMList(n),2})
        %             connectToM = RMList(n);
        %             break;
        %         end
        if (isfield(kcf, 'Product'))
            bonds = kcf.Product.Bonds;
        else
            bonds = kcf.bond2;
        end
        indexBond(i) = FindBondBetween2atoms(bonds, R, M(i));
        patternRMbonds(i) = bonds(indexBond(i), 4);
    else
        indexBond(i) = 0;
        patternRMbonds(i) = 0;
    end
end
end

function indexBond = FindBondBetween2atoms(bond, R, M)
% finding the index (column1) in bond where elements in col2 and 3 are R
% and M

arr1 = find(bond(:, 3) == R);
arr2 = find(bond(:, 2) == M);
indexBond = intersect(arr1,arr2);

if isempty(indexBond)
    arr1 = find(bond(:, 2) == R);
    arr2 = find(bond(:, 3) == M);
    indexBond = intersect(arr1,arr2);
end
end

function [cycle atoms bFound] = FindCycle(bonds, atom)
atoms = unique(bonds(:,2:3));
[cycle atoms bFound] = recursiveFn(bonds, [], [], atom);
end

function [cycle atoms bFound] = recursiveFn(bonds, cycle, atoms, atom)
bondIndx = FindNeighborsP(bonds, atom);
bFound = false;
for x = bondIndx.'
    if bonds(x, 3) == atom
        neighboringAtom = bonds(x,2);
    else
        neighboringAtom = bonds(x,3);
    end
    if ~isempty(atoms) && neighboringAtom==atoms(end), bFound=false; continue, end % all the leaves are disregarded here
    bFound = any(atom==atoms);
    if bFound
        cycle = [cycle x];
        atoms = [atoms atom];
        return
    else
        [cycleT atomsT bFoundT] = recursiveFn(bonds, [cycle x], [atoms atom], neighboringAtom);
        if bFoundT
            bFound=bFoundT;
            cycle = cycleT;
            atoms = atomsT;
            return
        end
    end
end
end

function neighbors = FindNeighborsP(bonds, atom)
col2 = find(bonds(:, 2) == atom);
col3 = find(bonds(:, 3) == atom);
neighbors = unique([col2; col3]);

end

function ConstructMolFile(inputStructure, filename)
nAtoms = size(inputStructure.Atoms, 1);
nBonds = size(inputStructure.Bonds, 1);
% constructing a mol file with atoms and bonds information
molFileText = [' \n \n \n ', num2str(nAtoms), ' ', num2str(nBonds),'  0  0  0  0  0  0  0  0999 V2000\n'];
for atom = 1:nAtoms
    limit = 9;
    nSpace = limit - length(num2str(inputStructure.Atoms{atom, 4}, '%.4f'));
    molFileText = [molFileText, repmat(' ', 1, nSpace)];
    molFileText = [molFileText, num2str(inputStructure.Atoms{atom, 4}, '%.4f')];
    limit = 10;
    nSpace = limit - length(num2str(inputStructure.Atoms{atom, 5}, '%.4f'));
    molFileText = [molFileText, repmat(' ', 1, nSpace)];
    molFileText = [molFileText, num2str(inputStructure.Atoms{atom, 5}, '%.4f'), '    0.0000  ', inputStructure.Atoms{atom, 3}, '  0 0 0 0 0 0 0 0 0 0 0 0\n'];
end
for bond = 1:nBonds
    if inputStructure.Bonds(bond,2)<10
        molFileText = [molFileText, '  ', num2str(inputStructure.Bonds(bond,2))];
    else
        molFileText = [molFileText, ' ', num2str(inputStructure.Bonds(bond,2))];
    end
    if inputStructure.Bonds(bond,3)<10
        molFileText = [molFileText, '  ' , num2str(inputStructure.Bonds(bond,3))];
    else
        molFileText = [molFileText, ' ' , num2str(inputStructure.Bonds(bond,3))];
    end
    molFileText = [molFileText, '  ', num2str(inputStructure.Bonds(bond,4)), '  0     0  0\n'];
end
molFileText = [molFileText, 'M  END\n'];
fileID = fopen(filename, 'w');
if fileID<0, disp('Error, opening the file'), return , end
fprintf(fileID, molFileText);
fclose(fileID);
end

function [x_rotated y_rotated] = RotateGroupofPoints(x, y, x_center, y_center, theta)
% create a matrix of these points, which will be useful in future calculations
v = [x;y];

% create a matrix which will be used later in calculations
center = repmat([x_center; y_center], 1, length(x));

% define a rotation matrix
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

% do the rotation...
s = v - center; % shift points in the plane so that the center of rotation is at the origin
so = R*s; % apply the rotation about the origin
vo = so + center; % shift again so the origin goes back to the desired center of rotation
% this can be done in one line as:
% vo = R*(v - center) + center

% pick out the vectors of rotated x- and y-data
x_rotated = vo(1,:);
y_rotated = vo(2,:);
end

function [inputList2, inputListNumbering2] = FindRMsAfterPhaseI(inputStructure)
% For the added atom in phase I (D), find the neighbors and put them in a row of inputList

% Bond: 1st&2nd cols-atom#, 3rd col is the number of bonds
% first col of inputList is main atom
k = 1; counter = 1;
inputList1(counter, k) = inputStructure.Atoms(inputStructure.DphaseI, 2);
inputListNumbering(counter, k) = inputStructure.DphaseI;
k = k+1;
% other cols are its neighbors
neighbors = FindNeighbors(inputStructure.Bonds, inputStructure.Atoms{inputStructure.DphaseI, 1});
for nbrL = 1:length(neighbors)
    inputList1(counter, k) = inputStructure.Atoms(neighbors(nbrL), 2);
    inputListNumbering(counter, k) = inputStructure.Atoms{neighbors(nbrL), 1};
    k = k + 1;
end
% while(k < 6)
%     inputList1{counter, k} = [];
%     inputListNumbering(counter, k) = 0;
%     k = k+1;
% end
counter = counter + 1;
if strcmp(inputList1(1, 1), 'O2x')
    for i = 1:length(neighbors)
        k = 1;
        inputList1(counter, k) = inputList1(1, i+1);
        inputListNumbering(counter, k) = neighbors(i);
        oxyNeighbors = FindNeighbors(inputStructure.Bonds, neighbors(i));
        k = k + 1;
        for nbrL = 1:length(oxyNeighbors)
            inputList1(counter, k) = inputStructure.Atoms(oxyNeighbors(nbrL), 2);
            inputListNumbering(counter, k) = inputStructure.Atoms{oxyNeighbors(nbrL), 1};
            k = k + 1;
        end
        counter = counter + 1;
    end
end
inputList = inputList1;
[inputList(:, 2:end), inputListNumbering(:, 2:end)] = SortRow(inputList1(:, 2:end), inputListNumbering(:, 2:end));
inputList2 = Augment(inputList);
inputListNumbering2 = Augment(inputListNumbering);
end

function [output1, output2] = SortRow(inputList, inputListNumbering)
output2 = [];
temp = inputList;
for i = 1:size(inputList, 1)
    index = ~cellfun('isempty', temp(i,:));
    [sorted, order] = sortrows(temp(i, index)');
    output1(i, 1:sum(index)) = sorted';
    if ~isempty(inputListNumbering)
        output2(i, 1:sum(index)) = inputListNumbering(i, order');
    end
end
end

function output = Augment(input)
% since the pattern has 4 neighbors, for any main with less than 4N, we
% augment the input with empty cells
[nMain nCol] = size(input);
nNeighbors = nCol - 1;
if isnumeric(input(1,1))
    output = zeros(nMain, 5);
else
    output{nMain, 5} = [];
end
if nNeighbors == 4
    output = input;
else
    output(:, 1:nCol) = input;
end
end

function [neighboringAtom neighborsPos] = FindNeighbors(bonds, atom)
neighboringAtom = [];neighborsPos=[];
if isempty(bonds)
    return
end
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