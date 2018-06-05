import string 
import os
# install pubchempy first and then import it to python as pcp
import pubchempy as pcp
import openbabel, pybel
import pandas as pd
from pandas import Series, DataFrame
    
from contextlib import contextmanager
import sys, os

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout
            
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
        with suppress_stdout():
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
 

def getCompoundsPubChemID(metName):
    compound = pcp.get_compounds(metName, 'name')
    for cmp in compound:
        return str(cmp.cid)
            
    
    
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
    molObj = pybel.readstring("mol", product)
    return molObj.formula
        
if __name__ == '__main__':
    # fileLocation = os.path.dirname(os.path.realpath(__file__))
    # molFilePath = os.path.join(".", "productsMol", "product_23.mol")
    # molFilePath = "product_10.mol"
    # CIDDetails = main_function(molFilePath)
    # print CIDs     
    # print CIDDetails
    
    # product = read_File(molFilePath)
    # molObj = pybel.readstring("mol", product)
    # dir(molObj.OBMol)
    # print molObj.formula
    # formula = getCompoundFormula(59721035)
    # print formula
    
    cid = getCompoundsPubChemID('CO2')
    print cid