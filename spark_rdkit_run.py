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



"""This file contains the code described in the journal article in
   Molecular Informatics: https://doi.org/10.1002/minf.201800082

The main guard largely contains code to parse commandline arguments and
setup the spark context. The experiment is run for a various amount of
instance sizes in order to compare their run times.

Run pyspark park_rdkit_run.py -h in order to get the number of
required arguments.

"""

import time
import os
import json
from datetime import datetime
import argparse

import numpy as np
from pyspark import SparkConf, SparkContext
from pyspark.sql import SQLContext
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

from utils import (rdd2smile, clean, comp_desc, desc_dict)


def run_experiment(filepath):
    """Run experiment for comparison on runtime for different amount of compounds as SMILES strings.

    The experiment is run for a various number block sizes.


    Parameters:
    ----------
    filepath : str path to parquet file
    """

    # list of descriptors from the RDKit Descriptors modules
    descriptors = list(np.array(Descriptors._descList)[:, 0])
    # the calculator module from RDKit
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptors)
    df = sql_context.read.parquet(filepath)
    sqlContext.read.parquet

    with open('logfile.json', 'w') as json_file:
        for num_lines in [1e5, 2e5, 5e5, 1e6, 1.5e6, 2e6]:
            smiles_rdd = df.select('can_smiles') \
                .limit(num_lines).rdd.repartition(10)
            # map the previously defined function to the SMILES RDD
            cleaned_rdd = smiles_rdd.map(clean)
            # convert rdd to data frame, remove None values, assign to cleaned_df
            cleaned_df = cleaned_rdd.map(
                lambda x: (x,)).toDF(['smiles']) \
                .dropna(thresh=1, subset=('smiles'))
            smiles_clean_rdd = cleaned_df.rdd.repartition(10)
            start_query = datetime.now()
            # create RDD of MOL by using the map module and
            # MolFromSmiles from RDKit, run cleaning again
            mol_rdd = smiles_clean_rdd.map(
                lambda x: Chem.MolFromSmiles(rdd2smile(x)))
            # assign results of the query for substructure 'CO' to sub_rdd
            # documentation from RDKit:
            # http://www.rdkit.org/Python_Docs/rdkit.Chem.rdchem.Mol-class.html#HasSubstructMatch
            sub_rdd = mol_rdd.map(
                lambda x: x.HasSubstructMatch(Chem.MolFromSmiles('CO')))
            count = sub_rdd.collect().count(True)
            end_query = datetime.now()
            desc_data_rdd = smiles_clean_rdd.map(
                lambda x: comp_desc(rdd2smile(x), calculator))
            descriptors_df = desc_data_rdd.map(
                lambda x: desc_dict(x)).toDF(descriptors)
            descriptors_df.write.parquet('descriptor_data_v'
                                         + str(time.time())
                                         + '.parquet')
            ex_dict = {
                'test time': str(datetime.now()),
                'molecules': str(num_lines),
                'totalTime': str(datetime.now() - start_time),
                'queryTime': str(end_query - start_query),
                'descTime': str(datetime.now() - end_query),
                'queryResultTrue': str(count)
            }

            json_file.dump(ex_dict)
            json_file.dump('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filepath', help='path to parquet file',
                        type=str)
    parser.add_argument('interpreter', help='interpreter to be used to'
                        + 'run the experiments.', type=str)
    args = parser.parse_args()
    path = args.filepath
    interpreter = args.interpreter

    os.environ['SPARK_MAJOR_VERSION'] = '1'
    os.environ['PYSPARK_PYTHON'] = interpreter

    conf = (SparkConf()
            .setAppName('app')
            .setMaster('yarn-client')
            .setAll([('spark.executor.memory', '80g'),
                     ('spark.executor.cores', '12'),
                     ('spark.cores.max', '12'),
                     ('spark.driver.memory', '80g')]))
    sc = SparkContext(conf=conf)
    sql_context = SQLContext(sc)

    start_time = datetime.now()
    run_experiment(path)
