#!/usr/bin/env python3
import numpy as np
class StrohSextic:
  '''
  Compute A, B, p in the Sextic Formalism of Stroh;
  Equations can be found in Chapter 5.1, 5.2, 5.3 [1]
  [1] TING, T. C. T. (1996). Anisotropic elasticity: Theory and applications. Oxford University Press. 
  '''
  def __init__(self,Cij=None,atol=1e-9):
    '''
    Cij should be in the current x1, x2, x3 coordinate system and specified as a 6x6 numpy array
    '''
    self.Cij=Cij
    self.atol=atol
    if self.Cij is not None: 
      self.sextic_equation()
  def calc(self,Cij=None,atol=1e-9):
    if Cij is not None:
      self.Cij=Cij
    self.tol=atol
    if self.Cij is not None:
      self.sextic_equation()
    else:
      raise Exception("Elastic modulus is not provided!!!")

  def sextic_equation(self):
    '''
    Define:
    A,B,p,S,H,L
    A,B are not unique;
    '''
    # Chapter 5.2 p137 [1]
    Q=np.array([[self.Cij[0,0],self.Cij[0,5],self.Cij[0,4]],
                [self.Cij[0,5],self.Cij[5,5],self.Cij[4,5]],
                [self.Cij[0,4],self.Cij[4,5],self.Cij[4,4]]])
    R=np.array([[self.Cij[0,5],self.Cij[0,1],self.Cij[0,3]],
                [self.Cij[5,5],self.Cij[1,5],self.Cij[3,5]],
                [self.Cij[4,5],self.Cij[1,4],self.Cij[3,4]]])
    T=np.array([[self.Cij[5,5],self.Cij[1,5],self.Cij[3,5]],
                [self.Cij[1,5],self.Cij[1,1],self.Cij[1,3]],
                [self.Cij[3,5],self.Cij[1,3],self.Cij[3,3]]])
    # Chapter 5.5 p144 [1]
    T_inv=np.linalg.inv(T)
    assert np.isclose(T_inv @ T,np.eye(3),rtol=0,atol=self.atol).all() and np.isclose(T @ T_inv,np.eye(3),rtol=0,atol=self.atol).all()
    N1=-T_inv @ R.T
    N2=T_inv
    N3=R @ T_inv @ R.T - Q
    N=np.block([[N1,N2],
                [N3,N1.T]])
    w,v=np.linalg.eig(N) # w: eigenvalues; v: eigenvectors;
    assert np.isclose(N@v,v@np.diag(w),rtol=0,atol=self.atol).all() # check whether Nv=wv
    assert not np.isclose(np.linalg.det(v),0,rtol=0,atol=self.atol) # eigenvectors should be linearly independent

    # order w and v; sort w in ascending order of positive imaginary part
    index=np.argsort(-np.imag(w))
    w=w[index]; v=v[:,index];
    index=np.array([2,1,0,3,4,5])
    w=w[index]; v=v[:,index];
    assert np.isclose(w[:3],np.conjugate(w[3:]),rtol=0,atol=self.atol).all() # w[:,3] should be the conjugate of w[3:]
    assert np.isclose(v[:,:3],np.conjugate(v[:,3:]),rtol=0,atol=self.atol).all() # v[:,:3] should be the conjugate of v[:,3:]

    A=v[:3,:3]; B=v[3:,:3]; p=w[:3]
    # normalize A, B such that B.T @ A + A.T @ B = I
    tmp=B.T @ A + A.T @ B # orthogonality relation
    assert np.isclose(tmp,np.diag(np.diag(tmp)),rtol=0,atol=self.atol).all() # tmp should be diagonal
    scale=np.diag(1/np.sqrt(np.diag(tmp),dtype='complex'))
    A=A @ scale
    B=B @ scale
    assert np.isclose(B.T @ np.conjugate(A) + A.T @ np.conjugate(B),np.zeros((3,3)),rtol=0,atol=self.atol).all() # orthogonality relation

    for i in range(3): # check (Q+p*(R+R.T)+p^2*T)a=0
      assert np.isclose((Q+p[i]*(R+R.T)+p[i]**2*T)@A[:,i],np.zeros(3),rtol=0,atol=self.atol).all()

    self.A=A; self.B=B; self.p=p
    # A dimension: Pa^-0.5; B dimension: Pa^0.5; p dimension: 1
    # Chapter 5.5 p146[1]
    S=1j * (2*self.A @ self.B.T - np.eye(3)); assert np.isclose(np.imag(S),0,rtol=0,atol=self.atol).all()
    H=2j * self.A @ self.A.T; assert np.isclose(np.imag(H),0,rtol=0,atol=self.atol).all()
    L=-2j * self.B @ self.B.T; assert np.isclose(np.imag(L),0,rtol=0,atol=self.atol).all()
    # S, H, L are all real, which can be derived from closure relations
    self.S=np.real(S) # dimension: 1
    self.H=np.real(H) # dimension: Pa^-1
    self.L=np.real(L) # dimension: Pa

    assert np.all(np.linalg.eigvalsh(self.L) > 0)
  def read_from_file(self,fname='stroh_tensor.dat'):
    '''
    Define:
    A,B,p
    '''
    fobj=open(fname,'r')
    StrohSextic._skip_comment(fobj)
    A_r=StrohSextic._read_matrix(fobj); A_i=StrohSextic._read_matrix(fobj);
    B_r=StrohSextic._read_matrix(fobj); B_i=StrohSextic._read_matrix(fobj);
    p_r=StrohSextic._read_matrix(fobj); p_i=StrohSextic._read_matrix(fobj);
    fobj.close()
    self.A=A_r+A_i*1j
    self.B=B_r+B_i*1j
    self.p=p_r+p_i*1j
    self.p=np.reshape(self.p,3)
    # Chapter 5.5 p146[1]
    S=1j * (2*self.A @ self.B.T - np.eye(3)); assert np.isclose(np.imag(S),0,rtol=0,atol=self.atol).all()
    H=2j * self.A @ self.A.T; assert np.isclose(np.imag(H),0,rtol=0,atol=self.atol).all()
    L=-2j * self.B @ self.B.T; assert np.isclose(np.imag(L),0,rtol=0,atol=self.atol).all()
    # S, H, L are all real, which can be derived from closure relations
    self.S=np.real(S) # dimension: 1
    self.H=np.real(H) # dimension: Pa^-1
    self.L=np.real(L) # dimension: Pa
  def _skip_comment(fobj):
    pos=fobj.tell()
    line=fobj.readline()
    while line != '':
      line=line.lstrip()
      if len(line)>0 and line[0] != '#':
        break
      pos=fobj.tell()
      line=fobj.readline()
    fobj.seek(pos)
  def _read_matrix(fobj):
    StrohSextic._skip_comment(fobj)
    line=fobj.readline()
    rc=line.strip().split()
    row=int(rc[0]); col=int(rc[1])
    ret=np.zeros((row,col))
    for i in range(row):
      StrohSextic._skip_comment(fobj)
      line=fobj.readline()
      line=line.strip().split()
      for j in range(col):
        ret[i,j]=float(line[j])
    return ret


  def to_file(self,fname='stroh_tensor.dat'):
    fobj=open(fname,'w')
    fmtstr='%.18e'
    fobj.write("# Stroh tensor A, B, p calculated using sextic formalism\n")
    fobj.write("# A real part:\n")
    fobj.write(str(self.A.shape[0])+' '+str(self.A.shape[1])+'\n')
    np.savetxt(fobj,np.real(self.A),fmtstr)
    fobj.write("# A imaginary part:\n")
    fobj.write(str(self.A.shape[0])+' '+str(self.A.shape[1])+'\n')
    np.savetxt(fobj,np.imag(self.A),fmtstr)
    fobj.write("# B real part:\n")
    fobj.write(str(self.B.shape[0])+' '+str(self.B.shape[1])+'\n')
    np.savetxt(fobj,np.real(self.B),fmtstr)
    fobj.write("# B imaginary part:\n")
    fobj.write(str(self.B.shape[0])+' '+str(self.B.shape[1])+'\n')
    np.savetxt(fobj,np.imag(self.B),fmtstr)
    fobj.write("# p real part:\n")
    fobj.write(str(self.p.shape[0])+' 1\n')
    np.savetxt(fobj,np.real(self.p),fmtstr)
    fobj.write("# p imaginary part:\n")
    fobj.write(str(self.p.shape[0])+' 1\n')
    np.savetxt(fobj,np.imag(self.p),fmtstr)
    fobj.close()
  def __eq__(self,other):
    '''
    Check equality of A,B,p; A, B are non-unique for a particular problem; reverse the sign of column vector doesn't affect physics
    Check equality of S,H,L
    Caution:
    where p[:] are not all unique the comparison algorithm may not work since the order of column vectors in A and B becomes undetermined
    '''
    atol=min(self.atol,other.atol)
    t1=np.all(np.isclose(self.p,other.p,rtol=0,atol=atol))
    for i in range(3):
      tmp1=np.all(np.isclose(self.A[:,i],other.A[:,i],rtol=0,atol=atol)) and \
           np.all(np.isclose(self.B[:,i],other.B[:,i],rtol=0,atol=atol))
      tmp2=np.all(np.isclose(self.A[:,i],-other.A[:,i],rtol=0,atol=atol)) and \
           np.all(np.isclose(self.B[:,i],-other.B[:,i],rtol=0,atol=atol))
      t1=t1 and (tmp1 or tmp2)
    t2=np.all(np.isclose(self.S,other.S,rtol=0,atol=atol)) and \
       np.all(np.isclose(self.H,other.H,rtol=0,atol=atol)) and \
       np.all(np.isclose(self.L,other.L,rtol=0,atol=atol))
    return t1 and t2


def test():
  from elastic_modulus import cij_hex
  from elastic_modulus import cij_cubic
  from elastic_modulus import ElasticModulus as em
  # Cij=cij_hex(C11=64.2662,C33=70.9312,C12=25.4545,C13=20.2856,C44=18.0165)
  Cij=cij_cubic(279.2,148.8,93.0)
  s=StrohSextic(Cij)
  # C=em(Cij)
  # newC=C.rotate_2(np.array([[0,0,1],[0,-1,0],[1,0,0]]))
  # newCij=newC.Cij
  # s=StrohSextic(newCij)
  s.to_file('test_data/stroh_tensor.dat')
  print("Write")
  print(s.p)
  print(s.A)
  print(s.B)
  print("Read")
  s.read_from_file('test_data/stroh_tensor.dat')
  print(s.p)
  print(s.A)
  print(s.B)

def main():
  test()

if __name__ == "__main__":
  main()
