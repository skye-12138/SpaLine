'''
This script is for ST data nucleus detected

'''
#!/usr/bin/env python
# conding:utf-8 -*-
import argparse
import stlearn as st
import scanpy as sc
import numpy as np
from numpy import random,mat
from pathlib import Path
import pandas as pd
from scipy import io,sparse
import os

def feature_tiling(inDir, outDir, sample):
    print (sample, "start SME normalize")

    #read data
    data = st.Read10X(path = inDir)
    data.var_names_make_unique()
    data.layers['raw_count'] = data.X
    #tile data
    TILE_PATH = Path(os.path.join(outDir,'{0}_tile'.format(sample)))
    TILE_PATH.mkdir(parents = True, exist_ok = True)

    #tile morphology
    st.pp.tiling(data,TILE_PATH, crop_size = 40)
    st.pp.extract_feature(data)

    return data




