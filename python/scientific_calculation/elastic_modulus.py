#!/usr/bin/env python3
import numpy as np
import re
import os
class ElasticModulus:
  '''
  Read in the elastic tensor in the 6x6 contracted notation
  References:
  [1] TING, T. C. T. (1996). Anisotropic elasticity: Theory and applications. Oxford University Press. 
  [2] Hirth, J. P. and Lothe J. (1982). Theory of dislocations.
  [3] Phillip L. Gould and Yuan Feng (2018) Introduction to Linear Elasticity
  '''
  def __init__(self,Cij,axes=np.eye(3,dtype='float64')):
    self.Cij=Cij
    self.axes=axes.astype('float64')
    ElasticModulus._normalize_axes(self.axes)
    ElasticModulus._check_axes(self.axes)
  def rotate_2(self,new_axes,tol=1e-10):
    '''
    Rotate Cij according to new_Cab = Kai Kaj Cij
    Formulas can be found in Chapter 2.8 p54 [1]
    new_axes: 3x3 numpy array; specified as [x',y',z']
    tol: set to 0 when Cij < tol
    '''
    new_axes=new_axes.astype('float64')
    ElasticModulus._normalize_axes(new_axes)
    ElasticModulus._check_axes(new_axes)
    R=np.matmul(new_axes,self.axes.T)
    K1=R**2
    K2=np.array([[R[0,1]*R[0,2],R[0,2]*R[0,0],R[0,0]*R[0,1]],
                 [R[1,1]*R[1,2],R[1,2]*R[1,0],R[1,0]*R[1,1]],
                 [R[2,1]*R[2,2],R[2,2]*R[2,0],R[2,0]*R[2,1]]])
    K3=np.array([[R[1,0]*R[2,0],R[1,1]*R[2,1],R[1,2]*R[2,2]],
                 [R[2,0]*R[0,0],R[2,1]*R[0,1],R[2,2]*R[0,2]],
                 [R[0,0]*R[1,0],R[0,1]*R[1,1],R[0,2]*R[1,2]]])
    K4=np.array([[R[1,1]*R[2,2]+R[1,2]*R[2,1],R[1,2]*R[2,0]+R[1,0]*R[2,2],R[1,0]*R[2,1]+R[1,1]*R[2,0]],
                 [R[2,1]*R[0,2]+R[2,2]*R[0,1],R[2,2]*R[0,0]+R[2,0]*R[0,2],R[2,0]*R[0,1]+R[2,1]*R[0,0]],
                 [R[0,1]*R[1,2]+R[0,2]*R[1,1],R[0,2]*R[1,0]+R[0,0]*R[1,2],R[0,0]*R[1,1]+R[0,1]*R[1,0]]])
    K=np.block([[K1,2*K2],[K3,K4]])
    new_Cij=np.einsum('ai,bj,ij->ab',K,K,self.Cij)
    new_Cij[abs(new_Cij)<tol]=0
    return ElasticModulus(new_Cij,new_axes)
  def rotate_4(self,new_axes,tol=1e-10):
    '''
    Rotate Cij by transforming it first to Cijkl, then doing the rotation, and then transforming it back to Cij
    new_axes: 3x3 numpy array; specified as [x',y',z']
    tol: set to 0 when Cij < tol
    '''
    new_axes=new_axes.astype('float64')
    ElasticModulus._normalize_axes(new_axes)
    ElasticModulus._check_axes(new_axes)
    R=np.matmul(new_axes,self.axes.T) # R[i,j]=new_axes[i,:] * axes[j,:]
    Cijkl=ElasticModulus.cij_to_cijkl(self.Cij)
    new_Cijkl=np.einsum('ai,bj,ck,dl,ijkl->abcd',R,R,R,R,Cijkl)
    new_Cij=ElasticModulus.cijkl_to_cij(new_Cijkl)
    new_Cij[abs(new_Cij)<tol]=0
    return ElasticModulus(new_Cij,new_axes)
  def set_axes(self,new_axes):
    self.axes=new_axes
  def vogit_ave(self):
    '''
    Formulas can be fould in Chapter 13.2 p.424 [2] and Chapter 4.3 p.85 [3]
    '''
    ciijj=self.Cij[0:3,0:3].sum()
    cijij=np.einsum('ii',self.Cij)+np.diagonal(self.Cij[3:,3:]).sum()
    lame=1/15*(2*ciijj-cijij)
    mu=1/30*(3*cijij-ciijj)
    E=mu*(3*lame+2*mu)/(lame+mu)
    nu=lame/(2*(lame+mu))
    K=E/(3*(1-2*nu))
    return (lame,mu,E,nu)
  def _vogit_ave(self):
    '''
    Formulas can be fould in Chapter 13.2 p.424 [2] and Chapter 4.3 p.85 [3]
    '''
    Cijkl=ElasticModulus.cij_to_cijkl(self.Cij)
    lame=1/15*(2*np.einsum('iijj',Cijkl)-np.einsum('ijij',Cijkl))
    mu=1/30*(3*np.einsum('ijij',Cijkl)-np.einsum('iijj',Cijkl))
    E=mu*(3*lame+2*mu)/(lame+mu)
    nu=lame/(2*(lame+mu))
    K=E/(3*(1-2*nu))
    return (lame,mu,E,nu)
  def vogit_cubic_ave(self):
    cubic_sym, C_cubic=is_cubic(self.Cij)
    if not cubic_sym:
      raise Exception("Not cubic symmetry in Cij !!!")
    Cij=cij_cubic(*C_cubic)
    H=2*Cij[3,3]+Cij[0,1]-Cij[0,0]
    lame=Cij[0,1]-H/5
    mu=Cij[3,3]-H/5
    E=mu*(3*lame+2*mu)/(lame+mu)
    nu=lame/(2*(lame+mu))
    K=E/(3*(1-2*nu))
    return (lame,mu,E,nu)
  def reuss_ave(self):
    '''
    Formulas can be fould in Chapter 13.2 p.424 [2] and Chapter 4.3 p.85 [3]
    '''
    Sij=np.linalg.inv(self.Cij)
    Sij[0:3,3:]/=2
    Sij[3:,0:3]/=2
    Sij[3:,3:]/=4
    siijj=Sij[0:3,0:3].sum()
    sijij=np.einsum('ii',Sij)+np.diagonal(Sij[3:,3:]).sum()
    E=15/(2*sijij+siijj)
    nu=1/15*(sijij-2*siijj)*E
    lame=E*nu/((1+nu)*(1-2*nu))
    mu=E/(2*(1+nu))
    K=E/(3*(1-2*nu))
    return (lame,mu,E,nu)
  def reuss_cubic_ave(self):
    '''
    Reuss assuming Cij has cubic symmetry
    '''
    cubic_sym, C_cubic=is_cubic(self.Cij)
    if not cubic_sym:
      raise Exception("Not cubic symmetry in Cij !!!")
    Cij=cij_cubic(*C_cubic)
    Sij=np.linalg.inv(Cij)
    J=Sij[0,0]-Sij[0,1]-Sij[3,3]/2
    E=1/(Sij[0,0]-2/5*J)
    nu=-(Sij[0,1]+J/5)*E
    lame=E*nu/((1+nu)*(1-2*nu))
    mu=E/(2*(1+nu))
    K=E/(3*(1-2*nu))
    return (lame,mu,E,nu)
  def vogit_reuss_hill_ave(self):
    '''
    Average K and mu from vogit and reuss; Match with experimental polycrystalline better
    From ref: 10.1088/0370-1298/65/5/307
    Converion formula in wikipedia: https://en.wikipedia.org/wiki/Elastic_modulus
    '''
    vogit=self.vogit_ave()
    reuss=self.reuss_ave()
    K_func=lambda E,nu : E/(3*(1-2*nu))
    K_V=K_func(vogit[2],vogit[3])
    K_R=K_func(reuss[2],reuss[3])

    K=(K_V+K_R)/2
    mu=(vogit[1]+reuss[1])/2

    E=9*K*mu/(3*K+mu)
    lame=K-2*mu/3
    nu=(3*K-2*mu)/2/(3*K+mu)
    return (lame,mu,E,nu)

  def stability(self):
    '''
    Test whether Cij is positive definite
    '''
    w,v=np.linalg.eigh(self.Cij)
    stable= np.all(w>0)
    unstable_direction=v[:,w<=0]
    return (stable,unstable_direction)
  def cij_to_cijkl(Cij):
    '''
    Transform the contracted notation to full notation
    Relevant formulas can be found in chapter 2.3 p36 [1]
    '''
    dict={0:[0,0],1:[1,1],2:[2,2],3:[1,2],4:[0,2],5:[0,1]}
    Cijkl=np.zeros((3,3,3,3))
    for i in range(6):
      for j in range(6):
        Cijkl[dict[i][0],dict[i][1],dict[j][0],dict[j][1]]=Cij[i,j]
        Cijkl[dict[i][0],dict[i][1],dict[j][1],dict[j][0]]=Cij[i,j]
        Cijkl[dict[i][1],dict[i][0],dict[j][0],dict[j][1]]=Cij[i,j]
        Cijkl[dict[i][1],dict[i][0],dict[j][1],dict[j][0]]=Cij[i,j]
    return Cijkl
  def cijkl_to_cij(Cijkl):
    dict={0:[0,0],1:[1,1],2:[2,2],3:[1,2],4:[0,2],5:[0,1]}
    Cij=np.zeros((6,6))
    for i in range(6):
      for j in range(6):
        Cij[i,j]=(Cijkl[dict[i][0],dict[i][1],dict[j][0],dict[j][1]]+\
                  Cijkl[dict[i][0],dict[i][1],dict[j][1],dict[j][0]]+\
                  Cijkl[dict[i][1],dict[i][0],dict[j][0],dict[j][1]]+\
                  Cijkl[dict[i][1],dict[i][0],dict[j][1],dict[j][0]])/4
    return (Cij+Cij.T)/2
  def _check_axes(axes,tol=1e-10):
    '''
    Check whether axes are orthogonal and whether it is a right hand coordinate system
    axes: 3x3 numpy array; [x,y,z]
    tol: x is considered orthogonal to y if np.dot(x,y)<tol
    '''
    if abs(np.dot(axes[0,:],axes[1,:])) > tol or abs(np.dot(axes[0,:],axes[2,:])) > tol or abs(np.dot(axes[1,:],axes[2,:])) > tol:
      raise Exception("Axes are not orthogonal to each other!!!")
    if np.linalg.det(axes) < 0:
      raise Exception("Not a right hand coordinate system!!!")
  def _normalize_axes(axes):
    '''
    normalize axes
    '''
    for i in range(3):
      axes[i,:]=axes[i,:]/np.linalg.norm(axes[i,:])

def read_cij_from_pot(fname,tol=1e-10):
  '''
  Read Cij from lammps calculation results
  Return Cij as a 6x6 numpy array
  '''
  if not os.path.exists(fname):
    return
  Cij=np.zeros((6,6))
  regex=re.compile(r'C([123456])([123456])\s+equal\s+(\S+)')
  with open(fname,'r') as f:
    for line in f:
      mo=regex.search(line)
      if mo:
        Cij[int(mo[1])-1,int(mo[2])-1]=float(mo[3])
        Cij[int(mo[2])-1,int(mo[1])-1]=float(mo[3])
  Cij[abs(Cij)<tol]=0
  return Cij
def cij_to_ci(Cij):
  '''
  Cij: (6,6) np.ndarray
  Return:
  Ci: (21,) np.ndarray
  Mapping:
  C1 , C7 , C8 , C9 , C10, C11
  C7 , C2 , C12, C13, C14, C15
  C8 , C12, C3 , C16, C17, C18
  C9 , C13, C16, C4 , C19, C20
  C10, C14, C17, C19, C5 , C21
  C11, C15, C18, C20, C21, C6 
  '''
  Ci=np.zeros(21)
  for i in range(6):
    Ci[i]=Cij[i,i]
  n=6
  for i in range(6):
      for j in range(i+1,6):
          Ci[n]=(Cij[i,j]+Cij[j,i])/2
          n+=1
  return Ci
def ci_to_cij(Ci):
  '''
  Ci: (21,) np.ndarray
  Return:
  Cij: (6,6) np.ndarray
  Mapping:
  C1 , C7 , C8 , C9 , C10, C11
  C7 , C2 , C12, C13, C14, C15
  C8 , C12, C3 , C16, C17, C18
  C9 , C13, C16, C4 , C19, C20
  C10, C14, C17, C19, C5 , C21
  C11, C15, C18, C20, C21, C6 
  '''
  Cij=np.zeros((6,6))
  for i in range(6):
      Cij[i,i]=Ci[i]
  n=6
  for i in range(6):
      for j in range(i+1,6):
          Cij[i,j]=Ci[n]
          n+=1
  Cij = np.where(Cij,Cij,Cij.T)
  return Cij

def is_cubic(Cij,tol_abs=1):
  '''
  Check whether Cij has cubic symmetry
  Return:
  [True|False, Cij_cubic]
  '''
  cubic_sym=True

  C11_idx=np.full((6,6),False)
  C11_idx[0,0]=C11_idx[1,1]=C11_idx[2,2]=True
  C11_all=Cij[C11_idx]

  C12_idx=np.full((6,6),False)
  C12_idx[0,1]=C12_idx[0,2]=C12_idx[1,2]=C12_idx[1,0]=C12_idx[2,0]=C12_idx[2,1]=True
  C12_all=Cij[C12_idx]

  C44_idx=np.full((6,6),False)
  C44_idx[3,3]=C44_idx[4,4]=C44_idx[5,5]=True
  C44_all=Cij[C44_idx]

  for C in [C11_all,C12_all,C44_all]:
    if (C.max()-C.min()) > tol_abs:
      cubic_sym=False

  zeros_idx=np.logical_not(C11_idx | C12_idx | C44_idx)
  zeros_all=Cij[zeros_idx]

  if np.any(np.abs(zeros_all) > tol_abs):
    cubic_sym=False

  C_cubic=[C11_all.mean(),C12_all.mean(),C44_all.mean()]

  return cubic_sym, C_cubic
def cij_cubic(C11,C12,C44):
  Cij=np.zeros((6,6))
  Cij[0,0]=Cij[1,1]=Cij[2,2]=C11
  Cij[0,1]=Cij[0,2]=Cij[1,2]=Cij[1,0]=Cij[2,0]=Cij[2,1]=C12
  Cij[3,3]=Cij[4,4]=Cij[5,5]=C44
  return Cij

def is_hex(Cij,tol_abs=1):
  '''
  Check whether Cij has hex symmetry
  Return:
  [True|False, Cij_hex]
  '''
  hex_sym=True

  C11_idx=np.full((6,6),False)
  C11_idx[0,0]=C11_idx[1,1]=True
  C11_all=Cij[C11_idx]

  C12_idx=np.full((6,6),False)
  C12_idx[0,1]=C12_idx[1,0]=True
  C12_all=Cij[C12_idx]

  C13_idx=np.full((6,6),False)
  C13_idx[0,2]=C13_idx[1,2]=C13_idx[2,0]=C13_idx[2,1]=True
  C13_all=Cij[C13_idx]

  C33_idx=np.full((6,6),False)
  C33_idx[2,2]=True
  C33_all=Cij[C33_idx]

  C44_idx=np.full((6,6),False)
  C44_idx[3,3]=C44_idx[4,4]=True
  C44_all=Cij[C44_idx]

  C66_idx=np.full((6,6),False)
  C66_idx[5,5]=True
  C66_all=np.array([Cij[5,5],(C11_all.mean()-C12_all.mean())/2])

  for C in [C11_all,C12_all,C13_all,C33_all,C44_all,C66_all]:
    if (C.max()-C.min()) > tol_abs:
      hex_sym=False

  zeros_idx=np.logical_not(C11_idx | C12_idx | C13_idx | C33_idx | C44_idx | C66_idx)
  zeros_all=Cij[zeros_idx]

  if np.any(np.abs(zeros_all) > tol_abs):
    hex_sym=False

  C_hex=[C11_all.mean(), C12_all.mean(), C13_all.mean(), C33_all.mean(), C44_all.mean()]
  return hex_sym, C_hex

def cij_hex(C11,C12,C13,C33,C44):
  Cij=np.zeros((6,6))
  Cij[0,0]=Cij[1,1]=C11
  Cij[0,1]=Cij[1,0]=C12
  Cij[0,2]=Cij[1,2]=Cij[2,0]=Cij[2,1]=C13
  Cij[2,2]=C33
  Cij[3,3]=Cij[4,4]=C44
  Cij[5,5]=(C11-C12)/2
  return Cij
def cij_iso(**kwargs):
  '''
  Provide any two of E, nu, lame, mu
  return Cij
  '''
  import sympy as sym
  lame_s,mu_s,E_s,nu_s=sym.symbols('lambda,mu,E,nu')
  eq1=sym.Eq(lame_s,E_s*nu_s/((1+nu_s)*(1-2*nu_s)))
  eq2=sym.Eq(mu_s,E_s/(2*(1+nu_s)))
  if 'E' in kwargs.keys() and 'nu' in kwargs.keys():
    sol=sym.solve([eq1,eq2],[lame_s,mu_s],dict=True)
    lame=sol[0][lame_s].subs({E_s:kwargs['E'],nu_s:kwargs['nu']})
    mu=sol[0][mu_s].subs({E_s:kwargs['E'],nu_s:kwargs['nu']})
  elif 'E' in kwargs.keys() and 'lame' in kwargs.keys():
    # this case is special because it has two solutions
    sol=sym.solve([eq1,eq2],[mu_s,nu_s],dict=True)
    lame=kwargs['lame']
    mu0=sol[0][mu_s].subs({E_s:kwargs['E'],lame_s:kwargs['lame']})
    mu1=sol[1][mu_s].subs({E_s:kwargs['E'],lame_s:kwargs['lame']})
    mu=mu1 if mu1 > mu0 else mu0
  elif 'E' in kwargs.keys() and 'mu' in kwargs.keys():
    sol=sym.solve([eq1,eq2],[lame_s,nu_s],dict=True)
    lame=sol[0][lame_s].subs({E_s:kwargs['E'],mu_s:kwargs['mu']})
    mu=kwargs['mu']
  elif 'nu' in kwargs.keys() and 'lame' in kwargs.keys():
    sol=sym.solve([eq1,eq2],[E_s,mu_s],dict=True)
    lame=kwargs['lame']
    mu=sol[0][mu_s].subs({nu_s:kwargs['nu'],lame_s:kwargs['lame']})
  elif 'nu' in kwargs.keys() and 'mu' in kwargs.keys():
    sol=sym.solve([eq1,eq2],[lame_s,E_s],dict=True)
    lame=sol[0][lame_s].subs({nu_s:kwargs['nu'],mu_s:kwargs['mu']})
    mu=kwargs['mu']
  elif 'lame' in kwargs.keys() and 'mu' in kwargs.keys():
    lame=kwargs['lame']
    mu=kwargs['mu']
  else:
    raise Exception("Enter any two of E, nu, lame, mu!!!")
  lame=sym.N(lame); mu=sym.N(mu)
  Cij=np.zeros((6,6))
  Cij[0,0]=Cij[1,1]=Cij[2,2]=lame+2*mu
  Cij[0,1]=Cij[0,2]=Cij[1,2]=Cij[1,0]=Cij[2,0]=Cij[2,1]=lame
  Cij[3,3]=Cij[4,4]=Cij[5,5]=mu
  return Cij

def test():
  import atomman as am
  import atomman.unitconvert as uc
  from sympy import Matrix, init_printing, pprint
  import time
  init_printing()
  # C11=uc.set_in_units(94.86,'GPa'); C12=uc.set_in_units(110.97,'GPa'); C44=uc.set_in_units(52.86,'GPa');
  C11=uc.set_in_units(256.47,'GPa'); C12=uc.set_in_units(152.01,'GPa'); C44=uc.set_in_units(144.55,'GPa');

  # testing for rotation
  x=np.array([1,1,1]); z=np.array([1,-1,0]); y=np.cross(z,x)
  axes=np.array([x,y,z])
  C_am=am.ElasticConstants(C11=C11,C12=C12,C44=C44)
  C_me=ElasticModulus(uc.get_in_units(C_am.Cij,'GPa'))
  print(C_me.vogit_ave())
  print(C_me._vogit_ave())
  print(C_me.vogit_cubic_ave())
  print(C_me.reuss_ave())
  print(C_me.reuss_cubic_ave())
  print(C_me.vogit_reuss_hill_ave())

  print("Current Program (rotate_2):")
  start=time.time()
  newC=C_me.rotate_2(axes)
  end=time.time()
  print("Runing time: "+str((end-start)*1e3)+'ms')
  Cij_me_2=newC.Cij
  # pprint(Matrix(Cij_me_2))

  print("Current Program (rotate_4):")
  start=time.time()
  newC=C_me.rotate_4(axes)
  print(newC.vogit_ave())
  print(newC.reuss_ave())
  print(C_me.vogit_reuss_hill_ave())
  end=time.time()
  print("Runing time: "+str((end-start)*1e3)+'ms')
  Cij_me_4=newC.Cij
  # pprint(Matrix(Cij_me_4))

  print("Result from atomman:")
  start=time.time()
  newC=C_am.transform(axes)
  end=time.time()
  print("Runing time: "+str((end-start)*1e3)+'ms')
  Cij_am=uc.get_in_units(newC.Cij,'GPa')
  # pprint(Matrix(Cij_am))
  assert np.isclose(Cij_me_2,Cij_me_4,rtol=0,atol=1e-10).all()
  assert np.isclose(Cij_me_2,Cij_am,rtol=0,atol=1e-10).all()

  # testing for cij_iso
  vogit=C_me.vogit_ave()
  lame=vogit[0]; mu=vogit[1]; E=vogit[2]; nu=vogit[3]
  Cijiso1=cij_iso(lame=lame,mu=mu)
  Cijiso2=cij_iso(E=E,nu=nu)
  Cijiso3=cij_iso(E=E,lame=lame)
  Cijiso4=cij_iso(E=E,mu=mu)
  Cijiso5=cij_iso(nu=nu,mu=mu)
  Cijiso6=cij_iso(nu=nu,lame=lame)
  assert np.isclose(Cijiso1,Cijiso2,rtol=0,atol=1e-10).all()
  assert np.isclose(Cijiso2,Cijiso3,rtol=0,atol=1e-10).all()
  assert np.isclose(Cijiso3,Cijiso4,rtol=0,atol=1e-10).all()
  assert np.isclose(Cijiso4,Cijiso5,rtol=0,atol=1e-10).all()
  assert np.isclose(Cijiso5,Cijiso6,rtol=0,atol=1e-10).all()

  print(C_me.stability())

def main():
  test()

if __name__ == "__main__":
  main()
