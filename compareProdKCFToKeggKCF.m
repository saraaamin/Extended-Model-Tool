
function prodCmpID = compareProdKCFToKeggKCF(prodKCFFile)
KEGGKcfFolder = '.\KEGGData\compound\compound\kcf\';
KEGGKcfFolderFiles = dir(KEGGKcfFolder);
       
prodCmpID = [];
index = 1;

for keggCompFileIdx = 3:length(KEGGKcfFolderFiles)
    currentKeggKCFFile = KEGGKcfFolderFiles(keggCompFileIdx).name;
    kkeggCmpKCFFile =  fileread([KEGGKcfFolder, currentKeggKCFFile]); 
    isKCFEqual = compare_two_kcfFiles(kkeggCmpKCFFile, prodKCFFile);
    
    %return all KEGG products that match the product file passed to the function 
    if isKCFEqual
        prodCmpID(index,1) = str2num(currentKeggKCFFile(2:6));
        index = index+1;
    end
end
            
end
