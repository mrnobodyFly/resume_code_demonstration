#!/usr/bin/env python3
import argparse
import re
import numpy as np
import atomman as am

from pathlib import Path
import sys

sys.path.append(str(Path.home() / 'template/python_pkg'))
from modified_atomman.Strain import Strain

class NonlinearDisplacement:
  '''
  Calculate polar coordinate (r,theta) and neighbor list using ref_fname
  Calculate excessive nonlinear displacement and strain field by comparing kload_fname and perfect_kload_fname
  Strain is computed by calculating the stretching and rotating of vectors connecting to neighbors
  neighbor list is computed 
  '''
  def __init__(self,kload_fname,perfect_kload_fname, xc, yc,
               ref_fname='dump.crack',cutoff=3):
    '''
    kload_fname: str; relaxed crack configuration
    perfect_kload_fname: str; crack configuration from elastic solution
    xc,yc: float; crack tip location
    ref_fname: str; unloaded crack
    cutoff: used in strain calculation; should be between first nearest neighbor and second nearest neighbor
    '''
    # load system
    self.undeformed_sys=am.load('atom_dump',ref_fname)
    self.kload_sys=am.load('atom_dump',kload_fname)
    self.perfect_kload_sys=am.load('atom_dump',perfect_kload_fname)
    # define variables
    self.xc=xc; self.yc=yc; self.cutoff=cutoff

    self.output_sys=self.kload_sys
    self.clear_properties()
  
  @property
  def r(self):
    if self.__r is None:
      self.__R=self.compute_polar()
    return self.__r
  @property
  def theta(self):
    if self.__theta is None:
      self.__R=self.compute_polar()
    return self.__theta
  @property
  def R(self):
    if self.__R is None:
      self.__R=self.compute_polar()
    return self.__R

  @property
  def strain(self):
    if self.__strain is None:
      self.__strain=self.compute_strain()
    return self.__strain
  @property
  def e_xx(self):
    if self.__e_xx is None:
      strain=self.strain
      self.__e_xx=strain[:,0,0]
    return self.__e_xx
  @property
  def e_yy(self):
    if self.__e_yy is None:
      strain=self.strain
      self.__e_yy=strain[:,1,1]
    return self.__e_yy
  @property
  def e_xy(self):
    if self.__e_xy is None:
      strain=self.strain
      self.__e_xy=strain[:,0,1]
    return self.__e_xy

  @property
  def strain_polar(self):
    if self.__strain_polar is None:
      self.__strain_polar=self.compute_polar_strain()
    return self.__strain_polar
  @property
  def e_tr(self):
    if self.__e_tr is None:
      strain_polar=self.strain_polar
      self.__e_tr=strain_polar[:,1,0]
    return self.__e_tr
  @property
  def e_tt(self):
    if self.__e_tt is None:
      strain_polar=self.strain_polar
      self.__e_tt=strain_polar[:,1,1]
    return self.__e_tt
  @property
  def e_tz(self):
    if self.__e_tz is None:
      strain_polar=self.strain_polar
      self.__e_tz=strain_polar[:,1,2]
    return self.__e_tz

  @property
  def disp(self):
    if self.__disp is None:
      self.__disp=self.compute_disp()
    return self.__disp
  @property
  def disp_polar(self):
    if self.__disp_polar is None:
      self.__disp_polar=self.compute_polar_disp()
    return self.__disp_polar

  def compute_polar(self):
    # r, theta
    new_pos=self.undeformed_sys.atoms.prop('pos')-np.array([self.xc,self.yc,0]);
    theta=np.arctan2(new_pos[:,1],new_pos[:,0]); r=np.sqrt(new_pos[:,0]**2+new_pos[:,1]**2);
    # rotation matrix
    cos_theta=new_pos[:,0]/r; sin_theta=new_pos[:,1]/r;
    R=np.zeros([new_pos.shape[0],3,3]);
    R[:,0,0]=R[:,1,1]=cos_theta; R[:,0,1]=sin_theta; R[:,1,0]=-sin_theta; R[:,2,2]=1
    self.__r=r; self.__theta=theta/np.pi*180
    return R

  def compute_strain(self):
    neighbors=am.NeighborList(system=self.undeformed_sys,cutoff=self.cutoff)
    # standard atomman strain
    # strain=am.defect.Strain(system=self.kload_sys,neighbors=neighbors,basesystem=self.perfect_kload_sys,baseneighbors=neighbors)
    # customized strain calculation
    strain=Strain(system=self.kload_sys,basesystem=self.perfect_kload_sys,neighbors=neighbors)
    return strain.strain
  def compute_polar_strain(self):
    R=self.R; strain=self.strain
    strain_polar=np.einsum('nai,nbj,nij->nab',R,R,strain)
    return strain_polar

  def compute_disp(self):
    disp=self.kload_sys.atoms.prop('pos')-self.perfect_kload_sys.atoms.prop('pos')
    return disp
  def compute_polar_disp(self):
    R=self.R; disp=self.disp
    disp_polar=np.einsum('nai,ni->na',R,disp)
    return disp_polar

  def save_to_system(self,properties=None):
    all_keys=[
      'r','theta','strain','strain_polar','disp','disp_polar','e_tr','e_tt','e_tz','e_xx','e_yy','e_xy'
    ]
    default_keys=[
      'r','theta','e_tr','e_tt','e_tz','e_xy','e_yy'
    ]
    if properties is None:
      properties=default_keys
    for p in properties:
      assert p in all_keys, 'unknown property ' + p
      self.output_sys.atoms.view[p]=getattr(self,p)

  def to_dump_file(self,ofname='dump.tmp',prop_name=None):
    compulsory_prop_name=[
      'atom_id','atype','pos'
    ]
    default_prop_name= compulsory_prop_name + [
      'v_group_label','r','theta','e_tr','e_tt','e_tz','e_xx','e_yy','e_xy'
    ]
    if prop_name is None:
      prop_name=default_prop_name
    else:
      prop_name=compulsory_prop_name+prop_name

    for p in prop_name:
      assert p in self.output_sys.atoms.prop(), 'unknown property ' + p

    self.output_sys.dump('atom_dump',f=ofname,prop_name=prop_name)

  def interpolate_plot(self,prop_name='e_tr',ofname=None,
                       lx=60,ly=60,title=None):
    if ofname is None:
      ofname=prop_name+'.pdf'
    if title is None:
      title=prop_name
    settings={}
    settings['title']=title
    settings['xlim']=(self.xc-lx/2,self.xc+lx/2); settings['ylim']=(self.yc-ly/2,self.yc+ly/2);
    settings['plotxaxis']=[1,0,0]; settings['plotyaxis']=[0,1,0]
    settings['cmap']='bwr'
    fig=am.plot.interpolate_contour(system=self.output_sys,prop_name=prop_name,**settings)[2]
    fig.savefig(ofname)

  def clear_properties(self):
    self.__r=None
    self.__theta=None
    self.__R=None
    self.__strain=None
    self.__e_xx=None
    self.__e_yy=None
    self.__e_xy=None
    self.__strain_polar=None
    self.__e_tr=None
    self.__e_tt=None
    self.__e_tz=None
    self.__disp=None
    self.__disp_polar=None

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
  regex=re.compile('dump.k_load_(-?\d+\.\d+)_(-?\d+\.\d+)_(-?\d+\.\d+)')
  mo=regex.search(fname)
  if mo:
    # kii=float(mo[1]); ki=float(mo[2]); kiii=float(mo[3])
    kii=mo[1]; ki=mo[2]; kiii=mo[3]
  else:
    raise Exception("Wrong fname format")
  return (kii,ki,kiii)

def test():
  nld=NonlinearDisplacement('dump.k_load_0.6500_0.8500_0.0000','dump.perfect_k_load_0.6500_0.8500_0.0000','dump.crack')
  nld.compute_polar_disp(302.781,299.581)
  nld.to_dump_file("dump.test")
    

def main():
  # ./nonlinear_displacement.py -f dump.k_load_0.6500_0.8500_0.0000 -r dump.crack -c ./data/crack_pos.mod
  parser=argparse.ArgumentParser(
    prog="python3 nonlinear_displacement.py",
    description="Compute nonlinear displacement",
    epilog="Written by LI Songwei")
  parser.add_argument('-f','--kload_fname',type=str,required=True,help="k load dump file name")
  parser.add_argument('-r','--ref_fname',type=str,default='dump.crack',help="reference dump file name")
  parser.add_argument('-c','--crack_pos_fname',type=str,default='./data/crack_pos.mod',help="crack position file")
  args=parser.parse_args()
  kii,ki,kiii=read_k_from_fname(args.kload_fname)
  perfect_kload_fname='dump.perfect_k_load_'+kii+'_'+ki+'_'+kiii
  xc,yc=crack_pos(args.crack_pos_fname)
  ofname='dump.nonlinear_disp_'+kii+'_'+ki+'_'+kiii

  nld=NonlinearDisplacement(kload_fname=args.kload_fname,perfect_kload_fname=perfect_kload_fname,xc=xc,yc=yc,
                            ref_fname=args.ref_fname,cutoff=4)
  nld.save_to_system(properties=['r','theta','e_tr','e_tt','e_xx','e_yy','e_xy'])
  nld.to_dump_file(ofname,prop_name=['r','theta','e_tr','e_tt','e_xx','e_yy','e_xy'])
  # nld.interpolate_plot(title=r'$\varepsilon_{\theta r}$')


if __name__ == "__main__":
  main()

# import ovito
# pipeline=ovito.pipeline.Pipeline()
# pipeline.source=ovito.pipeline.FileSource()
# pipeline.source.load('dump.k_load_0.6500_0.8500_0.0000')
# disp_mod=ovito.modifiers.CalculateDisplacementsModifier()
# disp_mod.reference=ovito.pipeline.FileSource()
# pipeline.modifiers.append(disp_mod)
# data=pipeline.compute()
# ovito.io.export_file(data,'dump.tmp','lammps/dump',columns=["Particle Identifier","Particle Type",'Position.X','Position.Y','Position.Z','v_group_label','Displacement.X','Displacement.Y','Displacement.Z','Displacement Magnitude'])
