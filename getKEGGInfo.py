# from bioservices import *
from Bio import KEGG
from Bio.KEGG import REST
import re

# Function to retreive KEGG ID of a metabolite
def get_metabolite_ID(metName):
	
    keggData = REST.kegg_find('compound', metName).read()
    keggData = "".join(keggData)
    keggData = keggData.lower()
    keggData = keggData.replace('cpd:', '')
    keggData = keggData.replace('\t', '; ')
    keggData = keggData.splitlines()
    keggDataList = [item.split('; ') for item in keggData]
    
    keggID = ''
    
    for cmpItem in keggDataList:
        if metName.lower() in cmpItem:
            keggID = cmpItem[0]
            break
    
    if keggID != '':
        keggID = int(keggID[1:])
    else: 
        keggID = 0
    
    return keggID    

# Function to map a list of genes to theis EC numbers
# for a specific orgnaism
def map_EC_to_org(org):

    # open the file from KEGG containing EC numbers details
    org = org + ':'
    ecNum_genes_dict = {}
    with  open('enzyme', 'r') as openedFile:
        contents = openedFile.read()
        ecEntries = contents.split('///')
        
        # extract all genes and associated EC numbers and save to a dictionary
        for ecEntry in ecEntries:
            ecLines = ecEntry.splitlines()
            if ecLines[0] == '':
                ecLines.pop(0)
            ecNum = re.search('(\d+.\d+.\d+.\d+)', ecLines[0]).group(0)
            genes = []

            for i in range(1, len(ecLines)):
                if ecLines[i].find(org) != -1:
                    genesItems = ecLines[i].split()
                    genesItems.pop(genesItems.index(org))
                    if 'GENES' in genesItems:
                        genesItems.pop(0)
                        
                    for gene in genesItems:
                        genes.append(gene[0:5])
            
            for gene in genes:
                
                if gene in ecNum_genes_dict.keys():
                    ecNumList = ecNum_genes_dict[gene]
                    ecNumList.append(ecNum)
                    ecNum_genes_dict[gene] = ecNumList 
                else:
                    ecNum_genes_dict[gene] = [ecNum]
    
    writeFile = open('iM1515_genes_EC.csv', 'wb')
    writeFile.truncate()
    
    # extract gene-EC numbers mapping for only genes in iML1515 model
    with open('iM1515 gene.txt', 'r') as openedFile:
        contents = openedFile.read().splitlines()
        for gene in contents:
            try:
                ecNumList = ecNum_genes_dict[gene]
                line = gene + ',' + ';'.join(ecNumList) + '\n'
            except KeyError:
                line = gene + ',' + '' + '\n'
                
            writeFile.write(line)
    writeFile.close()    


def get_EC_num(geneName):
	
    # retrieve gene data from KEGG
    keggData = REST.kegg_find('genes', geneName).read()
    keggData = "".join(keggData)
    keggData = keggData.lower()
    keggData = keggData.splitlines()
    
    # find which line 'eco' exists in the returned values to get enzyme name
    enzymeNameLine = ''
    for line in keggData:
        if line.find('eco:') != -1:
            enzymeNameLine = line
            break
    
    if enzymeNameLine == '':
        return ''
    else:    
        enzymeName = enzymeNameLine[enzymeNameLine.index('\t')+1: enzymeNameLine.index(';')]        
    
    # find enzyme name in KEGG and get ECNums associated to it
    keggData = REST.kegg_find('enzyme', enzymeName).read()
    keggData = "".join(keggData)
    keggData = keggData.lower()
    keggData = keggData.splitlines()
    
    ecNumList = []
    for line in keggData:
        try:
            ecNumList.append(line[line.index(':')+1 : line.index('\t')])
        except ValueError:
            return ''

    if ecNumList == []:
        return ''
    else:
        ecNumList = ','.join(ecNumList)
        return ecNumList
    
 
def get_inchiKey(keggID): 

    kegg_con = KEGG()
    kegg_entry = kegg_con.parse(kegg_con.get(keggID))
    chebi_con = ChEBI()

    try:
        chebi_entry = chebi_con.getCompleteEntity('CHEBI:' + kegg_entry['DBLINKS']['ChEBI'])
        print chebi_entry
        return chebi_entry.inchiKey
    except:
        return '0'

        
if __name__ == '__main__':
    metName = ['4-Phospho-D-erythronate', '(S)-3-Hydroxydodecanoyl-CoA', '(S)-3-Hydroxydecanoyl-CoA', '(S)-3-Hydroxyoctanoyl-CoA', 'D-Fructose', 'Acetaldehyde', '3-Oxodecanoyl-[acyl-carrier protein]', '3-Dehydroquinate', 'trans-Dodec-2-enoyl-CoA', 'trans-Oct-2-enoyl-CoA', '5-Formamido-1-(5-phospho-D-ribosyl)imidazole-4-carboxamide', 'Succinyl-CoA', '5,10-Methylenetetrahydrofolate', 'GTP', '1-dodecanoyl-sn-glycerol 3-phosphate', 'Hydrogen sulfide', 'L-Proline', '5-phosphoribosyl-5-carboxyaminoimidazole', 'Dodecanoyl-CoA (n-C12:0CoA)', 'Octanoyl-CoA (n-C8:0CoA)', 'N-Acetyl-L-glutamate', 'Sodium', 'Acetate', '(R)-3-hydroxy-cis-palm-9-eoyl-[acyl-carrier protein]', '(R)-3-hydroxy-cis-dodec-5-enoyl-[acyl-carrier protein]', 'R-3-hydroxypalmitoyl-[acyl-carrier protein]', 'chorismate', 'L-Arginine', 'trans-Dec-2-enoyl-CoA', '1-hexadecanoyl-sn-glycerol 3-phosphate', 'N-(5-Phospho-D-ribosyl)anthranilate', 'ADP-D-glycero-D-manno-heptose', 'D-Glutamate', 'ADP-L-glycero-D-manno-heptose', 'Isocitrate', 'D-Glucosamine 6-phosphate', '2-Hydroxy-3-oxopropanoate', 'UTP', 'Acetyl phosphate', 'Diphosphate', 'Ornithine', 'Succinate', '(R)-3-Hydroxyhexanoyl-[acyl-carrier protein]', '(3R)-3-Hydroxyacyl-[acyl-carrier protein]', '(R)-3-Hydroxydodecanoyl-[acyl-carrier protein]', '(R)-3-hydroxy-cis-myristol-7-eoyl-[acyl-carrier protein]', 'UDP-N-acetylmuramoyl-L-alanyl-D-glutamate', 'dTDP', 'Indole', 'IMP', '(S)-3-Hydroxyhexanoyl-CoA', 'Adenosine', 'Deoxyguanosine']
    for met in metName:
        keggID = get_metabolite_ID(met)
        print keggID
    
    # inchikey = get_inchiKey('C00001')
    # print inchikey
    
    # map_EC_to_org('ECO')