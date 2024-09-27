#!/usr/bin/env python3
import argparse
import numpy as np
import re
import ovito

from pathlib import Path
import os, sys
sys.path.append(str(Path.home() / 'template/python_pkg'))
from elasticity.crack import Crack
from elasticity.elastic_modulus import read_cij_from_pot

def k_load_crack(k0,dk,nstep,dumpfname,crack_pos,cij,orient):
  pipeline=ovito.io.import_file(dumpfname)
  crack=Crack(Cij=cij,orient=orient,crack_loc=crack_pos)
  for i in range(nstep):
    k=k0+i*dk
    data=pipeline.compute()
    disp=crack.displacement(k)
    stress=crack.stress(k,'cart')
    pos=data.particles.position[:,0:2]
    data.particles_.positions_=disp(pos)+data.particles.positions
    s=stress(pos)
    s_vogit=np.zeros((s.shape[0],6))
    s_vogit[:,0]=s[:,0,0]
    s_vogit[:,1]=s[:,1,1]
    s_vogit[:,2]=s[:,2,2]
    s_vogit[:,3]=s[:,1,2]
    s_vogit[:,4]=s[:,0,2]
    s_vogit[:,5]=s[:,0,1]
    data.particles_.create_property('s',data=s_vogit)
    ofname='dump.python_k_load_{kii:.4f}_{ki:.4f}_{kiii:.4f}'.format(kii=k[0],ki=k[1],kiii=k[2])
    ovito.io.export_file(data,ofname,'lammps/dump',columns=['Particle Identifier', 'Particle Type', 'Position.X','Position.Y','Position.Z', 'v_group_label','s.1','s.2','s.3','s.4','s.5','s.6'])

def read_crack_pos(fname='crack_pos.mod'):
  regex_x=re.compile('variable\s+Xc\s+equal\s(\S+)')
  regex_y=re.compile('variable\s+Yc\s+equal\s(\S+)')
  fobj=open(fname,'r')
  Xc=None; Yc=None
  for line in fobj:
    if Xc is None:
      mo=regex_x.search(line)
      if mo:
        Xc=float(mo[1])
    elif Yc is None:
      mo=regex_y.search(line)
      if mo:
        Yc=float(mo[1])
    else:
      break
  fobj.close()
  return np.array([Xc,Yc])

def test():
  crack_pos=read_crack_pos()
  k=np.array([0.33,0.4,0])
  dk=np.array([0.01,0,0])
  nstep=2
  dumpfname='dump.crack'
  cij=read_cij_from_pot('elastic_c_Ni_fcc_static_0K.mod')
  orient=np.loadtxt('orientation.dat')
  k_load_crack(k,dk,nstep,dumpfname,crack_pos,cij,orient)

def main():
   # ./k_load_crack.py -c 5 5 -n 2 -f dump.crack -o crack_pos.mod -t sharp -p 0 0 1
  parser=argparse.ArgumentParser(
    prog="python3 k_load_crack.py",
    description="Displace atoms according to elastic solution near crack tip",
    epilog="Written by LI Songwei")
  parser.add_argument('-k','--k0',type=float,nargs=3,required=True,help="Initial K")
  parser.add_argument('-d','--dk',type=float,nargs=3,required=True,help="dK loading step")
  parser.add_argument('-n','--nstep',type=int,help="Number of loading steps")
  parser.add_argument('-f','--dump_file',default='./dump.crack',type=str,help="Lammps dump file")
  parser.add_argument('-p','--crack_pos_file',default='./data/crack_pos.mod',type=str,help="crack tip position")
  parser.add_argument('-e','--cij_fname',type=str,help="Cij file name")
  parser.add_argument('-o','--orient_fname',type=str,help="Orientation file name")
  args=parser.parse_args()

  k0=np.array(args.k0)
  dk=np.array(args.dk)
  crack_pos=read_crack_pos(args.crack_pos_file)
  cij=read_cij_from_pot(args.cij_fname)
  orient=np.loadtxt(args.orient_fname)

  k_load_crack(k0,dk,args.nstep,args.dump_file,crack_pos,cij,orient)

  # kii=float(sys.argv[1])
  # ki=float(sys.argv[2])
  # kiii=float(sys.argv[3])
  # dkii=float(sys.argv[4])
  # dki=float(sys.argv[5])
  # dkiii=float(sys.argv[6])
  # nstep=int(sys.argv[7])
  # dumpcrack_file=sys.argv[8]
  # crack_pos_file=sys.argv[9]
  # cij_file=sys.argv[10]
  # orient_file=sys.argv[11]

  # k0=np.array([kii,ki,kiii])
  # dk=np.array([dkii,dki,dkiii])
  # crack_pos=read_crack_pos(crack_pos_file)
  # cij=read_cij_from_pot(cij_file)
  # orient=read_orientation(orient_file)

  # k_load_crack(k0,dk,nstep,dumpcrack_file,crack_pos,cij,orient)
if __name__ == "__main__":
  main()
