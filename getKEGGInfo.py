from Bio import KEGG
from Bio.KEGG import REST
from bioservices import *
import re

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
   
def map_EC_to_org(org):
    # open the file from KEGG containing EC numbers details
    org = org + ':'
    ecNum_genes_dict = {}
    with  open('enzyme', 'r') as openedFile:
        contents = openedFile.read()
        ecEntries = contents.split('///')
        
        for ecEntry in ecEntries:
            ecLines = ecEntry.splitlines()
            if ecLines[0] == '':
                ecLines.pop(0)
            print ('line', ecLines[0])
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
   
# def get_EC_num(geneName):
	
    # # retrieve gene data from KEGG
    # keggData = REST.kegg_find('genes', geneName).read()
    # keggData = "".join(keggData)
    # keggData = keggData.lower()
    # keggData = keggData.splitlines()
    
    # # find which line 'eco' exists in the returned values to get enzyme name
    # enzymeNameLine = ''
    # for line in keggData:
        # if line.find('eco:') != -1:
            # enzymeNameLine = line
            # break
    
    # if enzymeNameLine == '':
        # return ''
    # else:    
        # enzymeName = enzymeNameLine[enzymeNameLine.index('\t')+1: enzymeNameLine.index(';')]        
    
    # # find enzyme name in KEGG and get ECNums associated to it
    # keggData = REST.kegg_find('enzyme', enzymeName).read()
    # keggData = "".join(keggData)
    # keggData = keggData.lower()
    # keggData = keggData.splitlines()
    
    # ecNumList = []
    # for line in keggData:
        # try:
            # ecNumList.append(line[line.index(':')+1 : line.index('\t')])
        # except ValueError:
            # return ''

    # if ecNumList == []:
        # return ''
    # else:
        # ecNumList = ','.join(ecNumList)
        # return ecNumList
    
 
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
        
    # print chebi_entry.smiles
    # print chebi_entry.inchi
    
if __name__ == '__main__':
    # keggID = get_metabolite_ID('D-Glycerate 2-phosphate')
    # keggID = get_EC_num('s0001')
    # print keggID
    
    # inchikey = get_inchiKey('C01571')
    # print inchikey
    
    map_EC_to_org('ECO')