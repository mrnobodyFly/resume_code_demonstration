#!/usr/bin/env python3
import argparse
import numpy as np
from pathlib import Path

import os, sys
sys.path.append(str(Path.home() / 'template/python_pkg'))
from elasticity.elastic_modulus import read_cij_from_pot
from elasticity.stroh import StrohSextic
from elasticity.elastic_modulus import ElasticModulus

def main():
   # ./stroh_tensor.py -e Cij.mod -o orientation.dat -s stroh_tensor.dat
  parser=argparse.ArgumentParser(
    prog="python3 stroh_tensor.py",
    description="Calculate Stroh tensor",
    epilog="Written by LI Songwei")
  parser.add_argument('-e','--cij_fname',type=str,required=True,help="Cij file in lammps format")
  parser.add_argument('-o','--orient_fname',type=str,required=True,help="Orientation file")
  parser.add_argument('-s','--stroh_tensor_fname',default='./data/stroh_tensor.dat',type=str,help="Output Stroh tensor file")
  args=parser.parse_args()

  cij_fname=args.cij_fname
  orient_fname=args.orient_fname
  stroh_tensor_fname=args.stroh_tensor_fname

  new_axes=np.loadtxt(args.orient_fname)
  Cij=read_cij_from_pot(args.cij_fname)
  C=ElasticModulus(Cij)
  newC=C.rotate_2(new_axes)
  newCij=newC.Cij
  ss=StrohSextic(newCij)
  ss.to_file(args.stroh_tensor_fname)

  # # input argument:
  # # orient_fname cij_fname stroh_tensor_fname
  # orient_fname=sys.argv[1]
  # cij_fname=sys.argv[2]
  # stroh_tensor_fname=sys.argv[3]
  # # read files
  # # read orientation
  # new_axes=read_orientation(orient_fname)
  # # read cij
  # Cij=read_cij_from_pot(cij_fname)
  # # rotate cij according to orientation
  # C=ElasticModulus(Cij)
  # newC=C.rotate_2(new_axes)
  # newCij=newC.Cij
  # # calculate stroh tensor for rotated cij
  # ss=StrohSextic(newCij)
  # ss.to_file(stroh_tensor_fname)
if __name__ == "__main__":
  main()
