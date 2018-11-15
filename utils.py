'''
Author: Mario Lovric, Know-Center, Graz, Austria

MIT License
Copyright (c) 2018 Mario Lovric

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''


"""This file contains helper functions used in the experiment

"""

from rdkit import Chem


def rdd2smile(x):
    """converts SMILES as rdd's to SMILES as strings.

    Parameters:
    ----------
     x SMILES in form of RDD

    Returns:
    ----------
    SMILES in form string
    """
    return x[0].encode("ascii", "ignore")

def clean(y):
    """Uses the function rdd2smiles to convert RDD SMILES to SMILES
    string, if there is a value, otherwise it converts the SMILES
    string to a MOL file format assigned to the chx variable. The
    molecule in MOL format is subsequently converted back to an
    isomeric SMILES.

    Parameters
    ----------
    y RDD SMILES

    Returns
    ----------
    cleaned isomeric SMILES
    """
    if y is None:
        return None
    chx = Chem.MolFromSmiles(rdd2smile(y))
    if chx is None:
        return None
    return Chem.MolToSmiles(chx, isomericSmiles=True)

def comp_desc(smiles, calculator):
    """ Calculates descriptors for a given compunds as SMILES strings
    Parameters
    ----------
    smiles SMILES strings of a compound
    calculator the calculator module

    Returns
    ----------
    array of descriptor values for the compound
    """
    try:
        # assign converted structure from SMILES to MOL format
        mol = Chem.MolFromSmiles(smiles, sanitize=False) 
        # sanitize structure, module import from rdkit
        Chem.SanitizeMol(mol) 
        # assign result of descriptor calculation
        result = np.array(calculator.CalcDescriptors(mol)) 
         # RDKit returns -666 in case it is not able to return a valid value, please visit
         # https://github.com/rdkit/rdkit/search?q=-666&unscoped_q=-666
        if np.all(result == -666):
            # return an array of None values times the length of the descriptor list
            return np.array([np.nan]*len(descriptors)) 
        # if the array of descriptors is not finite,
        if np.all(np.isfinite(result)) == False: 
            # return an array of None values times the length of the descriptor list
            return np.array([np.nan]*len(descriptors)) 
        return result
    except:
        # if the functions throws an exception, return an array of None values
        return np.array([np.nan]*len(descriptors))

def desc_dict(x):
    """ Returns a dictionary where keys are descriptors names and values are calculated descriptors
    Parameters
    ----------
    x list of calculated descriptor values

    Returns
    ----------
    dictionary containing descriptors names and values
    """
    d = {}
    for i in range(len(x)):
        d[descriptors[i]] = float(x[i])
    return d
def _fng_mol(mol):
    """
    Parameters
    ----------
    smiles string, preg. canonical already
    Returns
    -------
    fingeprint mols if no typerror
    """

    try:
        return FingerprintMols.FingerprintMol(mol)
    except TypeError:
        return None
def _sim(l1,l2):
    """
    Parameters
    ----------
    l1 fingerprint mol
    l2 fingerprint mol
    Returns
    -------
    Similarity
    """
    try:
        return DataStructs.FingerprintSimilarity(l1,l2)
    except:
        return None
