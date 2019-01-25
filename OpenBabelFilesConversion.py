import string 
import os
# install pubchempy first and then import it to python as pcp
import pubchempy as pcp
import openbabel, pybel
import pandas as pd
from pandas import Series, DataFrame
            
def convert_mol_to_smiles(molFilePath, molFile):
    fileLocation = os.path.dirname(os.path.realpath(__file__))
    smilesFilePath = os.path.join(".", "productsSmiles", molFile)
    
    # conver mol file to smiles using opaenBable    
    # molFilePath_original = molFilePath
    
    smilesFilePath = smilesFilePath.replace('mol','smiles')
    # smilesFilePath = smilesFilePath.replace('productsMol','productsSmiles')
    # print molFilePath
    # print smilesFilePath
    
    os.system('babel -h -i mol '+ molFilePath + ' -o smi ' + smilesFilePath)
    return smilesFilePath    
    
def read_File(fileAddress):
    # read the smiles file
    f = open(fileAddress,'r')
    product_smiles = f.read().rstrip()
    return product_smiles
    
def search_by_SMILES(smiles):
    # search the data base
    CIDs = []
    CIDDetails = []
    try:
        results = pcp.get_compounds(smiles, 'smiles')
        for result in results: 
            # CIDs.append(result.cid)
            # CIDDetails.append(pcp.Compound.from_cid(result.cid).iupac_name)
            CIDs = result.cid
            CIDDetails = pcp.Compound.from_cid(result.cid).iupac_name
            break
    except:
        print 'serverError'
    
    # print CIDs
    # print CIDDetails
    if CIDs == None:
        details = ''
    else:
        details = str(CIDs) + ',' + str(CIDDetails)
    
    return details

def getCompoundFormula(cid):
    cidFormula = pcp.Compound.from_cid(cid).molecular_formula
    return cidFormula    
 
def getCompoundInchi(cid):
    cidInchi = pcp.Compound.from_cid(cid).inchi
    cidInchiKey = pcp.Compound.from_cid(cid).inchikey
    return cidInchi, cidInchiKey

def getCompoundSmiles(cid):
    cidSmiles = pcp.Compound.from_cid(cid).canonical_smiles
    cidIsomericsmiles = pcp.Compound.from_cid(cid).isomeric_smiles
    return cidSmiles, cidIsomericsmiles
    
def getCompoundsPubChemID(metName):
    compound = pcp.get_compounds(metName, 'name')
    for cmp in compound:
        return str(cmp.cid)
            
    
def getPubchemIDFromInchiKey(inchikey):
    compound = pcp.get_compounds(inchikey, 'inchikey')    
    for cmp in compound:
        return cmp.cid
        
    
def main_function(molFile):
    fileLocation = os.path.dirname(os.path.realpath(__file__))
    molFilePath = os.path.join(".", "productsMol", molFile)

    matches = {}
    path_to_smilesFile = convert_mol_to_smiles(molFilePath, molFile)
    product_smiles = read_File(path_to_smilesFile)
    CIDDetails = search_by_SMILES(product_smiles)
    
    return CIDDetails 

def convert_mol_to_formula(molFile):
    fileLocation = os.path.dirname(os.path.realpath(__file__))
    molFilePath = os.path.join(".", "productsMol", molFile)
    
    product = read_File(molFilePath)
    try:
        molObj = pybel.readstring("mol", product)
    except IOError:
        return ''
    return molObj.formula

# This function is to read a list of pubchem IDs from a file, get some details
# related to those IDs and save them to a file
def getPubChemIDsDetails():
    IDsFile = open('prodPubChemIDs.txt', 'r')
    pubChemIDsList = IDsFile.read().splitlines()
    
    writeFile = open('pubchemIDs_Details.txt', 'wb')
    writeFile.truncate()
    header = 'Pubchem ID, InChi, InChiKey, Canonical Smiles, Isomeric Smiles\n'
    writeFile.write(header)
    
    for pubchemID in pubChemIDsList:
        print pubchemID
        inchi, inchikey = getCompoundInchi(pubchemID)
        canonical_smiles, isomeric_smiles = getCompoundSmiles(pubchemID)
        line = pubchemID + ';' + inchi + ';' + inchikey + ';' + canonical_smiles + ';' + isomeric_smiles + '\n'
        writeFile.write(line)

    writeFile.close()
    
    
# This function is to read a list of pubchem IDs from a file, get some details
# related to those IDs and save them to a file
def getInchiKeyDetails():
    IDsFile = open('prodPubChemIDs.txt', 'r')
    inchikeyList = IDsFile.read().splitlines()
    
    writeFile = open('inchikey_Details.csv', 'wb')
    writeFile.truncate()
    header = 'Pubchem ID, InChiKey\n'
    writeFile.write(header)
    
    count = 0
    for inchikey in inchikeyList:
        count += 1 
        print count 
        pubchemID = getPubchemIDFromInchiKey(inchikey)
        line = str(pubchemID) + ',' + inchikey + '\n'
        writeFile.write(line)

    writeFile.close()    
    
if __name__ == '__main__':
    # fileLocation = os.path.dirname(os.path.realpath(__file__))
    # molFilePath = os.path.join(".", "productsMol", "product_23.mol")
    molFilePath = "product_11.mol"
    CIDDetails = main_function(molFilePath)
    # print CIDDetails
    
    # product = read_File(molFilePath)
    # molObj = pybel.readstring("mol", product)
    # dir(molObj.OBMol)
    # print molObj.formula
    # formula = getCompoundFormula(59721035)
    # print formula
    
    # cid = getCompoundsPubChemID('CO2')
    # print cid
    # inchi, inchikey = getCompoundInchi(76137728)
    # print inchi
    # print inchikey
    
    # cidSmiles, cidIsomericsmiles = getCompoundSmiles(76137728)
    # print cidSmiles
    # print cidIsomericsmiles
    getPubChemIDsDetails()
    # getPubchemIDFromInchiKey('RHVHUHHCHFFRKC-KODPNVEZSA-N')
    # getInchiKeyDetails()