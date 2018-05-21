'''
Author: Mario Lovric, Know-Center, Austria
'''
from datetime import datetime
startTime = datetime.now()
#imports
import numpy as np
import os
os.environ['SPARK_MAJOR_VERSION'] = "1"
os.environ['PYSPARK_PYTHON'] = "/usr/bin/python"

from pyspark import SparkConf, SparkContext
from pyspark.sql import SQLContext

conf=SparkConf().setAppName('app').setMaster('yarn-client').setAll([('spark.executor.memory', '80g'), ('spark.executor.cores', '12'), ('spark.cores.max', '12'), ('spark.driver.memory','80g')])
sc = SparkContext(conf=conf)
sqlContext =SQLContext(sc)
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import time


def rdd2smile(x):
    return x[0].encode("ascii", "ignore")
def clean(y):
    if y is None:
        return None
    chx=Chem.MolFromSmiles(rdd2smile(y))
    if chx is None:
        return None
    return Chem.MolToSmiles(chx, isomericSmiles=True)
def compDesc(smiles, calculator):
    try:
        mol=Chem.MolFromSmiles(smiles,sanitize=False)
        Chem.SanitizeMol(mol)
        result = np.array(calculator.CalcDescriptors(mol))
        if np.all(result==-666):
            return np.array([np.nan]*len(descriptors))
        if np.all(np.isfinite(result)) == False:
            return np.array([np.nan]*len(descriptors))
        return result
    except:
        return np.array([np.nan]*len(descriptors))
def desc_dict(x):
    d = {}
    for i in range(len(x)):
        d[descriptors[i]] = str(x[i])
    return d

#list of descirptors
descriptors = list(np.array(Descriptors._descList)[:,0])
#define descriptor calculator
calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptors)

#read smiles into DF
df=sqlContext.read.parquet('/user/mlovric/pubchem_2mil.parquet')

if os.path.exists('ex_time.txt'):
    os.remove('ex_time.txt')

for i in [100000,200000,500000,1000000, 1500000, 2000000, 2500000]:
    #create rdd, limit or not
    smiles_rdd=df.select('can_smiles').limit(i).rdd.repartition(10)
    # clean smiles
    cleaned_rdd = smiles_rdd.map(clean)
    # remove nans in df
    cleaned_df = cleaned_rdd.map(lambda x: (x,)).toDF(['smiles']).dropna(thresh=1, subset=('smiles'))
    # cleaned_df = cleaned_df.dropna(thresh=1,subset=('smiles'))
    # create clean rdd dropped nan
    smiles_clean_rdd = cleaned_df.rdd.repartition(10)
    #queries
    ######
    startQuery=datetime.now()
    ######
    mol_rdd = smiles_clean_rdd.map(lambda x: Chem.MolFromSmiles(rdd2smile(x)))
    sub_rdd = mol_rdd.map(lambda x: x.HasSubstructMatch(Chem.MolFromSmiles('CO')))
    Count=sub_rdd.collect().count(True)
    ######
    endQuery=datetime.now()
    ######

    #calculate descriptors
    desc_data_rdd=smiles_clean_rdd.map(lambda x: compDesc(rdd2smile(x),calculator))
    #create dataframe with descriptors
    descriptors_df = desc_data_rdd.map(lambda x: desc_dict(x)).toDF(descriptors)
    #write file
    descriptors_df.write.parquet('descriptor_data_v'+str(time.time())+'.parquet')
    ######
    savingTime=datetime.now()###
    ######
    F=open('ex_time.txt','a')
    exDict={'test time':str(datetime.now()),'molecules':str(i),\
            'totalTime':str(datetime.now() - startTime),\
            'queryTime':str(endQuery - startQuery),\
            'descTime':str(savingTime - endQuery), \
            'queryResultTrue':str(Count)}
    F.write(str(exDict) + '\n')
    F.close()


