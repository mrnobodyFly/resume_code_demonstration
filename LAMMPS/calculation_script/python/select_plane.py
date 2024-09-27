#!/usr/bin/env python3
import sys
import argparse
import numpy as np

import atomman as am
from atomic_model.crystal_plane import atomic_plane_position, projection, plane_intersection

class CrackGeometry:
  def __init__(self,ck_center,data_fname):
    '''
    ck_center: numpy vector of length 2; store x, y coordinate of crack tip
    data_fname: str; LAMMPS data file name
    '''
    self.ck_center=ck_center
    system=am.load("atom_data",data_fname)
    self.r=system.atoms.pos[:,0:2] # nx2 numpy array
    index=((self.r-self.ck_center)**2).sum(axis=1) < 50**2
    self.r=self.r[index,:]
  def crack_plane(self,ncrack):
    '''
    Calculate upper and lower crack plane position; compute exact Yc (ck_center[1])
    ncrack: int; number of atomic layers to be removed
    update Yc
    '''
    planes, count=atomic_plane_position(self.r,n=np.array([0,1]),tol=0.1)
    ck_plane=projection(self.ck_center,np.array([0,1]))
    ck_i=np.searchsorted(planes,ck_plane,side='right')-1 #planes[ck_i] lies just below ck_plane
    if ncrack%2 == 0:
      self.ck_center[1]=(planes[ck_i]+planes[ck_i+1])/2
      self.upper_crack_plane_pos=(planes[ck_i+ncrack//2]+planes[ck_i+ncrack//2+1])/2
      self.lower_crack_plane_pos=(planes[ck_i-ncrack//2+1]+planes[ck_i-ncrack//2])/2
    else:
      if (ck_plane-planes[ck_i]) > (planes[ck_i+1]-ck_plane):
        ck_i=ck_i+1
      self.ck_center[1]=planes[ck_i]
      self.upper_crack_plane_pos=(planes[ck_i+ncrack//2]+planes[ck_i+ncrack//2+1])/2
      self.lower_crack_plane_pos=(planes[ck_i-ncrack//2]+planes[ck_i-ncrack//2-1])/2
    return self.ck_center, self.upper_crack_plane_pos, self.lower_crack_plane_pos
  def tip_plane(self,tip_type,plane_normal=np.array([-1,0,0])):
    '''
    Calculate crack tip position; compute ck_center, tip_loc
    update ck_center
    '''
    if tip_type == "sharp" or tip_type == "blunt":
      self.plane_normal=plane_normal
      self.tip_loc=self._tip_plane_helper(self.plane_normal)
      self.ck_center[0]=self.tip_loc[0]
    elif tip_type == "plane":
      self.plane_normal=plane_normal
      self.tip_loc=self._tip_plane_helper(self.plane_normal)
      p1=np.array([np.array([0,1]),self.ck_center])
      p2=np.array([self.plane_normal[:2],self.tip_loc])
      self.ck_center=plane_intersection(p1,p2)
    elif tip_type == "angle":
      self.plane_normal=np.reshape(plane_normal,(2,3))
      tip1=self._tip_plane_helper(self.plane_normal[0,:])
      tip2=self._tip_plane_helper(self.plane_normal[1,:])
      p1=np.array([self.plane_normal[0,:2],tip1])
      p2=np.array([self.plane_normal[1,:2],tip2])
      self.tip_loc=plane_intersection(p1,p2)
      self.ck_center=self.tip_loc
    return self.tip_loc
  def _tip_plane_helper(self,plane_normal):
    '''
    Calculate intersection of tip plane with upper/lower crack plane
    update ck_center
    '''
    planes, count=atomic_plane_position(self.r,n=plane_normal[:2],tol=0.1)
    tip=projection(self.ck_center,plane_normal[:2])
    tip_i=np.searchsorted(planes,tip,side='right') # adjacent plane lie left to ck_center
    tip=(planes[tip_i]+planes[tip_i-1])/2
    p1=np.array([plane_normal[:2],tip*plane_normal[:2]/np.linalg.norm(plane_normal[:2])]) # tip plane
    if plane_normal[1] >= 0:
      p2=np.array([[0,1],[self.ck_center[0],self.upper_crack_plane_pos]]) # upper crack plane
    else:
      p2=np.array([[0,1],[self.ck_center[0],self.lower_crack_plane_pos]]) # lower crack plane
    ckp=np.array([[0,1],self.ck_center])
    self.ck_center=plane_intersection(p1,ckp) # update self.ck_center[0]
    return plane_intersection(p1,p2)

  def to_file(self,ofile):
    with open(ofile,'w') as f:
      f.write("variable Xc equal "+str(self.ck_center[0])+'\n')
      f.write("variable Yc equal "+str(self.ck_center[1])+'\n')
      f.write("variable upper_crack_plane_pos equal "+str(self.upper_crack_plane_pos)+'\n')
      f.write("variable lower_crack_plane_pos equal "+str(self.lower_crack_plane_pos)+'\n')
      f.write("variable crack_height equal "+str(self.upper_crack_plane_pos-self.lower_crack_plane_pos)+'\n')
      f.write("variable tip_x equal "+str(self.tip_loc[0])+'\n')
      f.write("variable tip_y equal "+str(self.tip_loc[1])+'\n')
      if self.plane_normal.size == 3:
        f.write("variable n11g equal "+str(self.plane_normal[0])+'\n')
        f.write("variable n12g equal "+str(self.plane_normal[1])+'\n')
        f.write("variable n13g equal "+str(self.plane_normal[2])+'\n')
      if self.plane_normal.size == 6:
        f.write("variable n11g equal "+str(self.plane_normal[0,0])+'\n')
        f.write("variable n12g equal "+str(self.plane_normal[0,1])+'\n')
        f.write("variable n13g equal "+str(self.plane_normal[0,2])+'\n')
        f.write("variable n21g equal "+str(self.plane_normal[1,0])+'\n')
        f.write("variable n22g equal "+str(self.plane_normal[1,1])+'\n')
        f.write("variable n23g equal "+str(self.plane_normal[1,2])+'\n')


def main():
   # ./select_plane.py -c 5 5 -n 2 -f dump.crack -o crack_pos.mod -t sharp -p 0 0 1
  parser=argparse.ArgumentParser(
    prog="python3 select_plane.py",
    description="Calculate accurate crack plane position",
    epilog="Written by LI Songwei")
  parser.add_argument('-c','--ck_center',type=float,nargs=2,required=True,help="Approximate x, y coordinate of crack tip")
  parser.add_argument('-n','--ncrack',type=int,required=True,help="Number of atom layers to be removed")
  parser.add_argument('-f','--fname',default='data.init_sc',type=str,help="Lammps data file")
  parser.add_argument('-o','--ofile',default='./data/crack_pos.mod',type=str,help="Output file")
  parser.add_argument('-t','--tip_type',default='sharp',type=str,choices=['sharp','blunt','plane','angle'],help="Crack tip geometry")
  parser.add_argument('-p','--plane_normal',type=float,default=[-1,0,0],nargs='*',help="Crack tip plane normal")
  args=parser.parse_args()

  ck_gm=CrackGeometry(np.array(args.ck_center),args.fname)
  ck_gm.crack_plane(args.ncrack)
  ck_gm.tip_plane(args.tip_type,np.array(args.plane_normal))
  ck_gm.to_file(args.ofile)

if __name__ == '__main__':
  main()
