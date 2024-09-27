import numpy as np
import re
from pathlib import Path
import os, sys
sys.path.append(str(Path.home() / 'template/python_pkg'))

from elasticity.elastic_modulus import read_cij_from_pot
from elasticity.stroh import StrohSextic
from elasticity.elastic_modulus import ElasticModulus

def read_orientation(fname):
  '''
  Read c1,c2,c3 from orientation file
  Return a numpy array [[c1],[c2],[c3]]
  '''
  if not os.path.exists(fname):
    raise Exception(fname+" does not exists!!!")
  orient=np.loadtxt(fname)
  return orient
  # orient=np.zeros((3,3))
  # regex=re.compile('variable\s+c([123])([123])\s+equal\s+[\'\"]?(-?\d*\.?\d*/?\d*\.?\d*)[\'\"]?')
  # with open(fname,'r') as f:
  #   for line in f:
  #     mo=regex.search(line)
  #     if mo:
  #       print(line)
  #       print(mo[3])
  #       orient[int(mo[1])-1,int(mo[2])-1]=float(mo[3])
  # print(orient)
  # return orient
        
