#!/usr/bin/env python3
import numpy as np
from elasticity.stroh import StrohSextic as sh
from elasticity.elastic_modulus import ElasticModulus as em
from scipy.optimize import fsolve

class CohesivePlaneModel:
  '''
  Predict the onset of instability for cohesive zone model
  '''
  def __init__(self,Cij,orient):
    '''
    Cij: 6x6 numpy array in standard orientation; units in GPa
    orient: 3x3 numpy array [x1,x2,x3]; x1: crack extension; x2: crack plane normal; x3: crack line;
    '''
    C=em(Cij); newC=C.rotate_2(orient)
    self.Cij=newC.Cij; self.orient=orient;
    sh_sex=sh(self.Cij); self.sh_sex=sh_sex;
    self.Lambda_inv=2*self.sh_sex.L; # Lambda=L^-1 /2; Lambda^-1=2*L
    self.Lambda=np.linalg.inv(self.Lambda_inv) # units: GPa^-1

  def energy_release_rate(self,k):
    '''
    Compute energy release rate according to G=K Lambda K for horizontal crack
    k: (3,) | (n,3) np.ndarray; units in Mpa m^0.5
    Return:
    G: float | (n,) np.ndarray; units in mJ/m^2
    '''
    if len(k.shape) == 1 and k.shape[0] == 3:
      return np.einsum('i,ij,j',k,self.Lambda,k) * 1e6  # units: mJ/m^2
    elif len(k.shape) == 2 and k.shape[1] == 3:
      return np.einsum('ni,ij,nj->n',k,self.Lambda,k) * 1e6  # units: mJ/m^2
    else:
      raise Exception("Wrong input array shape !!!")

  def set_cohesive_plane(self,theta=0,s=np.array([0,1,0]),phi=None):
    '''
    theta: scalar; angle between cohesive plane and crack plane; units: radians
    s: (3,) np.array; unit vector in slip direction
    phi: a local maximum of phi along s; units: mJ/m^2
    '''
    R=R_theta(theta); self.R=R
    Lambda_theta_inv=np.einsum("ai,bj,ij->ab",R,R,self.Lambda_inv) # units: GPa

    self.Lambda_theta=np.linalg.inv(Lambda_theta_inv) # units: GPa^-1
    self.s=s/np.linalg.norm(s)
    self.F=F_theta(self.sh_sex.B,self.sh_sex.p,theta) # units: 1
    self.p=np.einsum('i,ij,j',self.s,Lambda_theta_inv,self.s) # units: GPa; not the same p in Stroh formalism
    self.phi=phi # units: mJ/m^2

  def effective_k(self,k):
    '''
    Compute effective K according to Keff=F.K
    k: (3,) | (n,3) np.ndarray; units in Mpa m^0.5
    Return:
    k_eff: (3,) | (n,3) np.ndarray; units in Mpa m^0.5
    '''
    if len(k.shape) == 1 and k.shape[0] == 3:
      return self.F@k
    elif len(k.shape) == 2 and k.shape[1] == 3:
      return np.einsum('ai,ni->na',self.F,k)
    else:
      raise Exception("Wrong input array shape !!!")

  def resolved_keff(self,k):
    '''
    Compute resolved effective K according to rKeff=s.F.K
    k: (3,) | (n,3) np.ndarray; units in Mpa m^0.5
    Return:
    k_eff: float | (n,) np.ndarray; units in Mpa m^0.5
    '''
    if len(k.shape) == 1 and k.shape[0] == 3:
      return self.s@self.F@k
    elif len(k.shape) == 2 and k.shape[1] == 3:
      return np.einsum('a,ai,ni->n',self.s,self.F,k)
    else:
      raise Exception("Wrong input array shape !!!")

  def effective_energy(self,k):
    '''
    Calculate effective energy release rate according to (s.FK)^2/p,
    which should be equal to gamma_us in the case of dislocation emission
    k: (3,) | (n,3) np.ndarray; units in Mpa m^0.5
    Return:
    G_eff: float | (n,) np.ndarray; effective energy release rate; units in mJ/m^2
    '''
    if len(k.shape) == 1 and k.shape[0] == 3:
      return np.einsum('a,ai,i',self.s,self.F,k)**2/self.p * 1e6 # units: mJ/m^2
    elif len(k.shape) == 2 and k.shape[1] == 3:
      return np.einsum('a,ai,ni->n',self.s,self.F,k)**2/self.p * 1e6 # units: mJ/m^2
    else:
      raise Exception("Wrong input array shape !!!")

  def crack_open_energy_release_rate(self,k):
    '''
    Compute energy release rate according to (F.K) Lambda_theta (F.K) for kinked crack
    The expression is only exact for the case of theta=0 | F=I
    k: (3,) | (n,3) np.ndarray; units in Mpa m^0.5
    Return:
    G: float | (n,) np.ndarray; energy release rate; units mJ/m^2
    '''
    if len(k.shape) == 1 and k.shape[0] == 3:
      return np.einsum('ai,i,ab,bj,j',self.F,k,self.Lambda_theta,self.F,k) * 1e6  # units: mJ/m^2
    elif len(k.shape) == 2 and k.shape[1] == 3:
      return np.einsum('ai,ni,ab,bj,nj->n',self.F,k,self.Lambda_theta,self.F,k) * 1e6  # units: mJ/m^2
    else:
      raise Exception("Wrong input array shape !!!")

  def critical_resolved_keff(self):
    '''
    Compute critical s.FK
    Return
    kc_eff: resolved effect K along constrained path direction; units in MPa m^0.5
    '''
    return np.sqrt(self.p * self.phi * 1e-6) # units: MPa m^0.5

  def critical_k(self,mode,**kwargs):
    '''
    mode: str; func|ratio|value|linear; how constrain on K is specified
    '''
    if mode == "func":
      return self.critical_k_general(**kwargs)
    elif mode == "ratio":
      return self.critical_k_ratio(**kwargs)
    elif mode == "value":
      return self.critical_k_value(**kwargs)
    elif mode == "linear":
      return self.critical_k_linear(**kwargs)
    else:
      raise Exception("Unknown mode !!!")

  def critical_k_general(self,k_constrain_1=lambda k:k[0],k_constrain_2=lambda k:k[2],k0=[0,1,0]):
    '''
    s: (3,) np.array; slip direction
    phi: a local maximum of phi along s
    '''
    f_criterion=criterion_fk(self.s,self.F,self.p,self.phi)
    f = lambda k : np.array([f_criterion(k),k_constrain_1(k),k_constrain_2(k)])
    kc=fsolve(f,k0,xtol=1e-8) # initial guess is important
    assert np.isclose(f(kc),np.zeros(3),rtol=0,atol=1e-6).all()
    return kc

  def critical_k_ratio(self,k_ratio=np.array([0,1,0])):
    '''
    k_ratio: (3,) np.ndarray; a unit vector specifying direction of K
    '''
    k_mag=np.sqrt(self.p*self.phi)*1e-3 / np.einsum('i,ij,j',self.s,self.F,k_ratio)
    return k_mag * k_ratio

  def critical_k_value(self,k1=[0,0],k2=[2,0]):
    '''
    ki: (idx,val)
    '''
    linear_constrain=np.zeros((2,4))
    linear_constrain[0,k1[0]]=1; linear_constrain[0,3]=k1[1]
    linear_constrain[1,k2[0]]=1; linear_constrain[1,3]=k2[1]
    return self.critical_k_linear(linear_constrain)

  def critical_k_linear(self,linear_constrain):
    '''
    linear_constrain: (2,4) np.ndarray
    Ax=b
    '''
    A=np.zeros((3,3)); b=np.zeros(3)
    A[0,:]=self.F.T @ self.s; b[0]=np.sqrt(self.p*self.phi)*1e-3
    A[1:,:]=linear_constrain[:,:3]; b[1:]=linear_constrain[:,3]
    return np.linalg.inv(A) @ b

  def energy_opening_k(self,mode='ratio',**kwargs):
    '''
    criterion for crack opening using energy based method
    '''
    if mode == "ratio":
      return self.energy_opening_k_ratio(**kwargs)
    elif mode == "func":
      return self.energy_opening_k_general(**kwargs)
    else:
      raise Exception("Unknown mode !!!")

  def energy_opening_k_ratio(self,k_ratio=np.array([0,1,0])):
    '''
    k_ratio: (3,) np.ndarray; a unit vector specifying direction of K
    '''
    k_mag=np.sqrt(self.phi*1e-6/np.einsum('ai,i,ab,bj,j',self.F,k_ratio,self.Lambda_theta,self.F,k_ratio))
    return k_mag * k_ratio # units: MPa m^0.5

  def energy_opening_k_general(self,k_constrain_1=lambda k:k[0],k_constrain_2=lambda k:k[2],k0=[10,10,10]):
    Lambda_theta=np.einsum('ai,bj,ij->ab',self.R,self.R,np.linalg.inv(self.sh_sex.L)/2)
    f_criterion=lambda k:np.einsum('ai,i,ab,bj,j',self.F,k,Lambda_theta,self.F,k)*1e3-self.phi*1e-3
    f = lambda k : np.array([f_criterion(k),k_constrain_1(k),k_constrain_2(k)])
    kc=fsolve(f,k0,xtol=1e-8) # initial guess is important
    assert np.isclose(f(kc),np.zeros(3),rtol=0,atol=1e-6).all()
    return kc

def criterion_fk(s,F,p,phi):
  '''
  instability criterion for the cohesive plane model
  s: (3,) np.array; unit 1; constrained path direction
  F: 3x3 np.array; unit 1; transformation matrix for K^eff
  p: scalar; unit GPa; p=s^T Lambda^-1s
  phi: scalar; unit mJ/m^2; maximum cohesive energy along the constrained path
  k: (3,) np.array; stress intensity factor; units in MPa m^0.5
  '''
  f=lambda k:np.einsum('i,ij,j',s,F,k)**2/p * 1e3 - phi * 1e-3 # units: J/m^2
  return f

def F_theta(B,p,theta):
  '''
  Matrix for calculating effective SIF. K^eff=F.K
  B: as defined in Stroh formalism in GPa^0.5
  p: as defined in Stroh formalism in 1
  theta: cohesive plane plane tilt angle in radians
  '''
  s=np.sin(theta); c=np.cos(theta)
  M=np.array([[-s*c,c**2-1/2,0],
              [s**2,-s*c,0],
              [0,0,-s]])
  N=np.array([[c**2-1/2,s*c,0],
              [-s*c,c**2,0],
              [0,0,c]])
  B_inv=np.linalg.inv(B)
  tmp=np.sqrt(c+p*s) # Which root should be used ???; branch cut at angle=-pi; argument in (-pi/2,pi/2]
  X=-np.real(B @ np.diag(p/tmp) @ B_inv)
  Y=np.real(B @ np.diag(1/tmp) @ B_inv)
  F=M@X+N@Y
  return F
def R_theta(theta):
  '''
  Rotation matrix
  theta: float; units in radians
  '''
  R=np.array([[np.cos(theta),np.sin(theta),0],
              [-np.sin(theta),np.cos(theta),0],
              [0,0,1]])
  return R

def calculate_theta_phi(orient,slip):
  slip_plane=np.cross(orient[2,:],slip)
  r_dir=np.cross(slip_plane,orient[2,:])
  angle=lambda v1,v2:np.arccos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))
  phi=np.pi/2-angle(slip,orient[2,:]); theta=np.pi/2-angle(r_dir,orient[1,:])
  return (theta,phi)

def test_hcp():
  from elasticity.elastic_modulus import cij_hex
  Cij=cij_hex(C11=64.2662,C33=70.9312,C12=25.4545,C13=20.2856,C44=18.0165)
  print("Testing using Mg/wu_2014_msmse")

  k=1.62293549696233 # c/a=1.62293549696233
  # basal/pyramidal I (0001)/[1-210]
  param={"orient":np.array([[0,1,0],[0,0,1],[1,0,0]]),"gamma":568,"gamma_us":319,"slip":np.array([-0.5,np.sqrt(3)/2,k])}
  theta,phi=calculate_theta_phi(param['orient'],param['slip']); s=np.array([np.cos(phi),0,np.sin(phi)])
  cpm=CohesivePlaneModel(Cij,param['orient'])
  cpm.set_cohesive_plane(0,np.array([0,1,0]),2*param['gamma'])
  # kc=cpm.critical_k('value',k1=[0,0],k2=[2,0])
  # kc=cpm.critical_k_ratio()
  # kc=cpm.critical_k_general()
  # kc=cpm.critical_k_value(k1=[0,0],k2=[2,0])
  kc=cpm.energy_opening_k('func')
  # kc=cpm.energy_opening_k_ratio()
  # kc=cpm.energy_opening_k_general()
  cpm.set_cohesive_plane(theta,s,param['gamma_us'])
  ke=cpm.critical_k_ratio()
  # ke=cpm.critical_k_general()
  # ke=cpm.critical_k_value(k1=[0,0],k2=[2,0])
  kc[abs(kc) < 1e-8]=0; ke[abs(ke) < 1e-8]=0
  print("--------------------------------")
  print("basal/pyramidal I (0001)/[1-210]: ")
  print("Kc = "+str(kc))
  print("ref = 0.255")
  print("Ke = "+str(ke))
  print("ref = 0.351")

def main():
  test_hcp()

if __name__ == "__main__":
  main()
