#!/usr/bin/env python3
import argparse
import numpy as np
import re

import atomman as am

from pathlib import Path
import sys

sys.path.append(str(Path.home() / 'template/python_pkg'))

class AtommanContourPlot:
  def __init__(self,dump_fname,xc,yc):
    self.atomic_system=am.load("atom_dump",dump_fname)
    self.xc=xc; self.yc=yc
  def plot(self,prop_name='e_tr',ofname=None,lx=40,ly=40,title=None):
    assert prop_name in self.atomic_system.atoms.prop(), 'unknown property' + prop_name 
    if ofname is None:
      ofname=prop_name+'.pdf'
    if title is None:
      title=prop_name
    settings={}
    settings['title']=title
    settings['xlim']=(self.xc-lx/2,self.xc+lx/2); settings['ylim']=(self.yc-ly/2,self.yc+ly/2);
    settings['plotxaxis']=[1,0,0]; settings['plotyaxis']=[0,1,0]
    settings['cmap']='bwr'
    fig=am.plot.interpolate_contour(system=self.atomic_system,prop_name=prop_name,**settings)[2]
    fig.tight_layout()
    fig.savefig(ofname)

def crack_pos(fname='./data/crack_info.mod'):
  regex=re.compile('variable\s+(Xc|Yc)\s+equal\s+(\S+)')
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
  # regex=re.compile('dump.nonlinear_disp_(-?\d+\.\d+)_(-?\d+\.\d+)_(-?\d+\.\d+)')
  regex=re.compile('dump.k_load_(-?\d+\.\d+)_(-?\d+\.\d+)_(-?\d+\.\d+)')
  mo=regex.search(fname)
  if mo:
    # kii=float(mo[1]); ki=float(mo[2]); kiii=float(mo[3])
    kii=mo[1]; ki=mo[2]; kiii=mo[3]
  else:
    raise Exception("Wrong fname format")
  return (kii,ki,kiii)

def main():
  # ./atomman_interpolate_contour_plot.py -f dump.nonlinear_disp_0.6500_0.8500_0.0000 -p e_tr -c ./data/crack_pos.mod
  parser=argparse.ArgumentParser(
    prog="python3 atomman_interpolate_contour_plot.py",
    description="Plot plastic strain",
    epilog="Written by LI Songwei")
  parser.add_argument('-f','--dump_fname',type=str,required=True,help="dump file name")
  parser.add_argument('-p','--prop_name',type=str,default='e_tr',help="property name")
  parser.add_argument('-c','--crack_pos_fname',type=str,default='./data/crack_pos.mod',help="crack position file")
  args=parser.parse_args()
  kii,ki,kiii=read_k_from_fname(args.dump_fname)
  dump_fname='dump.nonlinear_disp_'+kii+'_'+ki+'_'+kiii
  xc,yc=crack_pos(args.crack_pos_fname)
  ofname=args.prop_name+'_'+kii+'_'+ki+'_'+kiii+'.pdf'

  amcp=AtommanContourPlot(dump_fname,xc,yc)
  amcp.plot(prop_name=args.prop_name,ofname=ofname)

if __name__ == "__main__":
  main()
