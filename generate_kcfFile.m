%  this is a function that converts a mol file to a KCF file using the kegg
%  api
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
