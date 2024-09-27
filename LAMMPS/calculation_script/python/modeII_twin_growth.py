#!/usr/bin/env python3
import numpy as np
import glob
import ovito
import re
import argparse

pipeline=ovito.pipeline.Pipeline()
pipeline.source=ovito.pipeline.FileSource()
cna_mod=ovito.modifiers.CommonNeighborAnalysisModifier(mode=ovito.modifiers.CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff)
pipeline.modifiers.append(cna_mod)

def crack_pos(fname='./data/crack_pos.mod'):
  regex=re.compile(r'variable\s+(Xc|Yc)\s+equal\s+(\S+)')
  with open(fname,'r') as f:
    for line in f:
      mo=regex.search(line)
      if mo:
        if mo[1]=='Xc':
          Xc=float(mo[2])
        elif mo[1]=='Yc':
          Yc=float(mo[2])
  return (Xc, Yc)
def read_k_from_fname(fname):
  regex=re.compile(r'dump.k_load_(-?\d+\.\d+)_(-?\d+\.\d+)_(-?\d+\.\d+)')
  mo=regex.search(fname)
  if mo:
    kii=float(mo[1]); ki=float(mo[2]); kiii=float(mo[3])
  else:
    raise Exception("Wrong fname format")
  return (kii,ki,kiii)

def get_dump_fname():
  fnames=glob.glob('dump.k_load_*_0.0000_0.0000')
  kiis=[]
  for f in fnames:
    k=read_k_from_fname(f)
    kiis.append(k[0])
  tmp=zip(kiis,fnames)
  fnames=[f for _,f in sorted(zip(kiis,fnames),key=lambda pair:pair[0])]
  kiis.sort()
  return kiis, fnames

def main():
  xc,yc=crack_pos()
  kiis,fnames=get_dump_fname()
  twin_length=[]
  pipeline.source.load(fnames)
  for frame in range(pipeline.source.num_frames):
    data=pipeline.compute(frame)
    msk=(data.particles['v_group_label']==2) & (data.particles.structure_types!=2) & (data.particles.positions[:,1]>200) & (data.particles.positions[:,1]<225)
    twin_tip=data.particles.positions[msk,0].max()
    twin_length.append(twin_tip-xc)
  ret=np.array([kiis,twin_length]).T
  np.savetxt('twin_growth.dat',ret)

if __name__ == "__main__":
  main()
