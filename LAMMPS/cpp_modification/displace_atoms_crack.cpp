// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "displace_atoms_crack.h"

#include "atom.h"
#include "atom_vec_body.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "irregular.h"
#include "lattice.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "random_park.h"
#include "variable.h"

#include <cmath>
#include <cstring>

#include <complex>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>

using namespace LAMMPS_NS;
using MathConst::DEG2RAD;
using MathConst::MY_2PI;

/*
  linear elastic displacement field for decoupled in-plane anti-plane anisotropic materials
  subject to K=[KII,KI,KIII] loading

  displace_atoms_crack group-ID KII KI KIII Xc Yc stroh_file_path
*/

namespace {
  void read_stroh_tensor(char fname[],Eigen::Matrix3cd &A,Eigen::Matrix3cd &B,Eigen::Vector3cd &p);
  void read_matrix(std::istream &inStream, Eigen::MatrixXd & m);
  void skip_comment(std::istream & inStream);
}

/* ---------------------------------------------------------------------- */

DisplaceAtomsCrack::DisplaceAtomsCrack(LAMMPS *_lmp) : Command(_lmp) {}

/* ---------------------------------------------------------------------- */

DisplaceAtomsCrack::~DisplaceAtomsCrack() {}

/* ---------------------------------------------------------------------- */

void DisplaceAtomsCrack::command(int narg, char **arg)
{
  int i;

  if (domain->box_exist == 0)
    error->all(FLERR,"Displace_atoms_crack command before simulation box is defined");
  if (narg < 6) error->all(FLERR,"Illegal displace_atoms_crack command");
  if (modify->nfix_restart_peratom)
    error->all(FLERR,"Cannot displace_atoms_crack after "
               "reading restart file with per-atom info");

  if (comm->me == 0) utils::logmesg(lmp,"Displacing atoms near crack tip ...\n");

  // group and style

  igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find displace_atoms_crack group ID");
  groupbit = group->bitmask[igroup];

  if (modify->check_rigid_group_overlap(groupbit))
    error->warning(FLERR,"Attempting to displace atoms near crack tip in rigid bodies");

  double KII = atof(arg[1]); double KI = atof(arg[2]); double KIII = atof(arg[3]);
  double Xc = atof(arg[4]); double Yc = atof(arg[5]);

  using std::complex_literals::operator""i;
  Eigen::Matrix3cd A, B, B_inv, P;
  Eigen::Vector3cd p;
  Eigen::Vector3d K, u;
  read_stroh_tensor(arg[6],A,B,p);

  B_inv=B.inverse();
  P << 1, 0, 0,
    0, 1, 0,
    0, 0, 1;
  K << KII, KI, KIII;


  double pi = 3.14159265358979323846 ;
  // Each processor loop over the nlocal atoms it owns:
  
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) { // if atoms belong to the group specified in the command

      // calc delx, dely, delz        <------------------
      double x1 = x[i][0]-Xc;
      double x2 = x[i][1]-Yc;

      P(0,0)=std::sqrt(x1+p(0,0)*x2);
      P(1,1)=std::sqrt(x1+p(1,0)*x2);
      P(2,2)=std::sqrt(x1+p(2,0)*x2);

      u=1e7*std::sqrt(2*1e-10/pi) * (A*P*B_inv).real() * K;

      x[i][0] += u(0,0);
      x[i][1] += u(1,0);
      x[i][2] += u(2,0);
    }
  }

  // move atoms back inside simulation box and to new processors
  // use remap() instead of pbc() in case atoms moved a long distance
  // use irregular() in case atoms moved a long distance

  imageint *image = atom->image;
  for (i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->reset_box();
  auto irregular = new Irregular(lmp);
  irregular->migrate_atoms(1);
  delete irregular;
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // check if any atoms were lost

  bigint natoms;
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms && comm->me == 0)
    error->warning(FLERR,"Lost atoms via displace_atoms: original {} "
                   "current {}",atom->natoms,natoms);
}

namespace {
  void read_stroh_tensor(char fname[],Eigen::Matrix3cd &A,Eigen::Matrix3cd &B,Eigen::Vector3cd &p){
    std::ifstream inStream;
    inStream.open(fname);
    if (inStream.fail()){
      std::cerr << "Can not open " << fname << " !!!\n";
      exit(1);
    }
    // using std::complex_literals::operator""i;
    std::complex<double> i(0.0,1.0);
    Eigen::MatrixXd m;
    read_matrix(inStream,m);
    A=m;
    read_matrix(inStream,m);
    // A=A+m*1i;
    A=A+m*i;
    read_matrix(inStream,m);
    B=m;
    read_matrix(inStream,m);
    // B=B+m*1i;
    B=B+m*i;
    read_matrix(inStream,m);
    p=m;
    read_matrix(inStream,m);
    // p=p+m*1i;
    p=p+m*i;
  }

  void read_matrix(std::istream &inStream, Eigen::MatrixXd & m){
    skip_comment(inStream);
    int nrow, ncol;
    inStream >> nrow >> ncol;
    inStream.ignore(10000,'\n');
    m.resize(nrow,ncol);
    for (int i=0;i<nrow;i++){
      skip_comment(inStream);
      for (int j=0;j<ncol;j++){
        inStream >> m(i,j);
      }
      inStream.ignore(10000,'\n');
    }
  }

  void skip_comment(std::istream & inStream){
    // ignore comment lines start with '#'
    char nextSymbol;
    while (inStream >> nextSymbol){
      if (nextSymbol == '#')
        inStream.ignore(10000,'\n');
      else{
        inStream.putback(nextSymbol);
        break;
      }
    }
  }
}
