#!/usr/bin/env python3
import numpy as np
import os
from pathlib import Path
import glob
import warnings
import re
import pickle
from scipy.interpolate import make_interp_spline
import atomman as am

class PathMaker:
    '''
    Create Spline Interpolator for NEB path.
    '''

    def __init__(self, pkfname: str | None = None):
        if pkfname is None:
            self.box = None # am.Box; atomman box instance
            self.pbc = None # boolean np.ndarray
            self.atype = None # int np.ndarray
            self.natoms = None # int
            self.path_interpolator = None # created by make_interp_spline
            self.path_fnames = None # list
            self.spline_k = None # int
            self.reaction_coord = None # np.ndarray
            self.fitting_quality = None # dictionary: {"fit_error":np.ndarray,"tangent_similarity":,"test_error":}
        else:
            self.load(pkfname)

    def fit_path_interpolator(self, fnames, ftype='dump', sort_key = None, path_dir = '.', \
                              rc = 'accumulate',
                              pbc=None, spline_k=5, check_fit_error=False):
        '''
        Arguments:
            fnames: glob pattern | list of explicit name
            ftype: str; data|dump
            sort_key: function for sorting fnames
            path_dir: folder containing the path files
            rc: str|list; list of float: should be ascending
            pbc: list of boolean
            spline_k: B-spline degree; the first k-1 derivatives are continuous
            check_fit_error: boolean
        Action / Side effect:
            self.path_fnames
            self.reaction_coord
            self.path_interpolator
            self.spline_k
            self.box
            self.pbc
            self.atype
            self.natoms
        Return:
        '''
        # process glob pattern and sort fnames
        fnames=self._process_fnames(fnames=fnames,sort_key=sort_key,path_dir=path_dir)
        self.path_fnames = [ os.path.abspath(fn) for fn in fnames ]
        nimages = len(fnames)
        if nimages < spline_k+1:
            raise Exception(
                f'At leat {spline_k+1} images is required for spline interpolation!!!')
        # load image configuration
        reaction_coord,image_config =\
            self._load_image_config(fnames=fnames,ftype=ftype,rc=rc,pbc=pbc)
        self.reaction_coord=reaction_coord
        # spline fit image configuration with respect to reaction coordinate
        self.path_interpolator = make_interp_spline(
            reaction_coord, image_config, k=spline_k)
        self.spline_k=spline_k
        # check fitting quality
        if check_fit_error:
            self._check_fit_stability(reaction_coord,image_config)

    def _process_fnames(self,fnames,sort_key=None,path_dir='./'):
        '''
        Process glob pattern and sort fnames
        '''
        # process glob pattern
        if isinstance(fnames, str):
            fnames = glob.glob(os.path.join(path_dir,fnames))
        else:
            for i in range(len(fnames)):
                fnames[i]=os.path.join(path_dir,fnames[i])
        # sort fnames
        if sort_key is None:
            # regex_pattern = r'(?:dump|data)\.\S*_(\d+\.?\d*)$'
            sort_key=lambda fn : float(re.compile(r'(?:dump|data)\.\S*\.(\d+)$').search(fn)[1])
        fnames.sort(key=sort_key)
        return fnames

    def _load_image_config(self,fnames,ftype,rc,pbc=None):
        '''
        Load image configurations / path
        Arguments:
            fnames: list; lammps dump/data files in ascending order
            ftype; str; dump|data
            rc: str|list; 'accumulate'|'direct'|list of float
            pbc: boolean list
        Action / Side effect:
            self.box
            self.pbc
            self.atype
            self.natoms
        Return:
            reaction_coord_accumulate_norm: np.ndarray; reaction coordinate used by lammps NEB calculation
            image_config: (nconfig,natoms,3) np.ndarray; unwrapped atom position
            reaction_coord_direction_norm: np.ndarray; reaction coordinate defined as 2-norm distance from the initial configuration
        '''
        init_system = am.load('atom_'+ftype, data=fnames[0])
        self.box = init_system.box
        self.atype = init_system.atoms.atype
        self.natoms = init_system.atoms.natoms
        if ftype == 'data':
            if pbc is None:
                warnings.warn(
                    'Boundary condition is not given for data file !!! Assume ppp boundary condition!!!')
            else:
                init_system.pbc = pbc
        self.pbc = init_system.pbc

        nimages = len(fnames)
        image_config = np.zeros((nimages, self.natoms, 3))
        reaction_coord_direct_norm = np.zeros(nimages)
        reaction_coord_accumulate_norm = np.zeros(nimages)  # used by lammps
        image_config[0, :, :] = init_system.atoms.prop('pos')
        inter_system1 = init_system
        for i in range(1, nimages):
            inter_system2 = am.load('atom_'+ftype, data=fnames[i])
            if ftype == 'data':
                inter_system2.pbc = self.pbc
            dX = am.displacement(inter_system1, inter_system2,
                                 box_reference='initial') # minimum image criterion is applied
            inter_system1 = inter_system2
            image_config[i, :, :] = image_config[i-1, :, :]+dX
            reaction_coord_accumulate_norm[i] = reaction_coord_accumulate_norm[i -
                                                                               1] + np.linalg.norm(dX)
            reaction_coord_direct_norm[i] = np.linalg.norm(
                image_config[i, :, :]-image_config[0, :, :])
  
        if np.any(np.diff(reaction_coord_direct_norm) < 0):
            wraning.warn("Direct reaction coodination does not increase monotonically!!!")
        reaction_coord_direct_norm = reaction_coord_direct_norm / \
            reaction_coord_direct_norm[-1]
        reaction_coord_accumulate_norm = reaction_coord_accumulate_norm / \
            reaction_coord_accumulate_norm[-1]

        if type(rc) is str:
            if rc == 'accumulate':
                reaction_coord=reaction_coord_accumulate_norm
            elif rc == 'direct':
                reaction_coord=reaction_coord_direct_norm
            else:
                raise Exception("Unknown rc type!!!")
        elif type(rc) is list and len(rc) == nimages:
            reaction_coord=np.array(rc)
        else:
            raise Exception("Unknown rc type!!!")
        return reaction_coord, image_config

    def _check_fit_stability(self, rd, image_config):
        '''
        Check quality of the fitting
            rd: np.ndarray; reaction coordinate
            image_config: np.ndarray; image configuration
        Side effect:
             self.fitting_quality
        '''
        n = rd.shape[0]
        fit_error = np.zeros(n)
        for i in range(n):
            dx=self.path_interpolator(rd[i]) - image_config[i,:,:]
            fit_error[i]=np.linalg.norm(dx)

        tangent_similarity = np.ones(n)
        for i in range(1, n-1):
            dx = image_config[i+1, :, :]-image_config[i-1, :, :]
            t = self.path_interpolator(rd[i], 1)
            tangent_similarity[i] = (dx*t).sum() / \
                (np.linalg.norm(dx)*np.linalg.norm(t))

        test_error = np.zeros(n)
        for i in range(1, n-1):
            idx = np.full(n, True)
            idx[i] = False
            tmp_interpolator = make_interp_spline(
                rd[idx], image_config[idx, :, :], k=self.spline_k)
            dx = tmp_interpolator(rd[i])-image_config[i, :, :]
            test_error[i] = np.linalg.norm(dx)
        self.fitting_quality={"fit_error":fit_error,"tangent_similarity":tangent_similarity,"test_error":test_error}

        print("Fit error (0 for perfect fitting; should be nearly 0):")
        print(fit_error)
        print("Tangent error (1 for perfect fitting; should be >0.85):")
        print(tangent_similarity)
        print("k-1 Test norm error (0 for perfect fitting; should be <0.1):")
        print(test_error)
        if fit_error.max() > 1e-9 or tangent_similarity.min() < 0.85 or test_error.max() > 0.1: 
            warnings.warn("Fitting quality is not good!!!")
        # return fit_error, tangent_similarity, test_error

    def make_path(self,r,odir='.',ofname=None,replicate=None,scale=None,oftype='data'):
        '''
        r: float; range: 0-1
        odir: str; output directory
        ofname: str;
        replicate: None|list;
        scale: None|(3,3)np.ndarray; 
        '''
        if r < self.reaction_coord.min() or r > self.reaction_coord.max():
            raise Exception('Reaction coordinate outside fitting range!!!')
        if ofname is None:
            ofname = f'{oftype}.path_{r:.3f}'
            r = float(f'{r:.3f}')
        ofname = os.path.join(odir, ofname)

        # scaled interpolation
        if scale is None:
            scale = np.eye(3)
        X = self.path_interpolator(r) @ scale.T
        box = am.Box(vects=self.box.vects @ scale.T, origin=self.box.origin @ scale.T)
        atoms = am.Atoms(atype=self.atype, pos=X)
        system = am.System(atoms=atoms, box=box,
                           pbc=self.pbc, scale=False)
        # replicate system
        if not replicate is None:
            system = system.supersize(a_size=replicate[0],b_size=replicate[1],c_size=replicate[2])

        # write to file
        with open(ofname, 'w') as f:
            if oftype == 'data':
                system.dump('atom_data', f=f, atom_style='atomic', units='metal')
            elif oftype == 'dump':
                # system.dump('atom_dump', f=f, lammps_units='metal',scale=False)
                system.dump('atom_dump', f=f, lammps_units='metal')


    def make_pafiv1_path(self, r, odir='.', ofname=None, replicate=None, scale=None):
        '''
        r: float; range: 0-1
        odir: str; output directory
        ofname: str;
        replicate: None|list;
        scale: None|(3,3)np.ndarray; 
        '''
        if r < 0 or r > 1:
            raise Exception('Reaction coordinate must be between 0 and 1!!!')
        if ofname is None:
            ofname = 'data.pafiv1_'+f'{r:.3f}'
            r = float(f'{r:.3f}')
        ofname = os.path.join(odir, ofname)

        # scaled interpolation
        if scale is None:
            scale = np.eye(3)
        X = self.path_interpolator(r) @ scale.T
        dX = self.path_interpolator(r, 1) @ scale.T
        ddX = self.path_interpolator(r, 2) @ scale.T
        dddX = self.path_interpolator(r, 3) @ scale.T
        box = am.Box(vects=self.box.vects @ scale.T, origin=self.box.origin @ scale.T)
        atoms = am.Atoms(atype=self.atype, pos=X)
        system = am.System(atoms=atoms, box=box,
                           pbc=self.pbc, scale=False)
        system.atoms.dX = dX
        system.atoms.ddX = ddX
        system.atoms.dddX = dddX 
        # replicate
        if not replicate is None:
            system = system.supersize(a_size=replicate[0],b_size=replicate[1],c_size=replicate[2])
        # write to file
        with open(ofname, 'w') as f:
            system.dump('atom_data', f=f, atom_style='atomic')
        natoms = system.natoms
        pafi_info = np.zeros((natoms,10))
        pafi_info[:,0] = np.arange(1, natoms+1) # id: 1 to N
        pafi_info[:, 1:4] = system.atoms.dX # dX0
        normsq = (pafi_info[:, 1:4]**2).sum()
        pafi_info[:, 4:7] = system.atoms.ddX / normsq   # ddX0/|dX0|^2
        pafi_info[:, 7:10] = system.atoms.dddX / normsq # dddX0/|dX0|^2
        with open(ofname, 'a') as f:
            fmt = '%u'+' %.13g'*9
            np.savetxt(f, pafi_info, fmt=fmt, delimiter=' ',
                       header='\nPafiPath\n', comments='')

    def make_pafi_path(self, r, odir='.', ofname=None, replicate=None, scale=None):
        '''
        r: float; range: 0-1
        ofname: str;
        replicate: None|list;
        scale: None|(3,3)np.ndarray; 
        '''
        if r < 0 or r > 1:
            raise Exception('Reaction coordinate must be between 0 and 1!!!')
        if ofname is None:
            ofname = 'data.pafi_'+f'{r:.3f}'
            r = float(f'{r:.3f}')
        ofname = os.path.join(odir, ofname)

        # scaled interpolation
        if scale is None:
            scale = np.eye(3)
        X = self.path_interpolator(r) @ scale.T
        dX = self.path_interpolator(r, 1) @ scale.T
        ddX = self.path_interpolator(r, 2) @ scale.T
        box = am.Box(vects=self.box.vects @ scale.T, origin=self.box.origin @ scale.T)

        atoms = am.Atoms(atype=self.atype, pos=X)
        system = am.System(atoms=atoms, box=box,
                           pbc=self.pbc, scale=False)
        system.atoms.dX = dX
        system.atoms.ddX = ddX

        if not replicate is None:
            system = system.supersize(a_size=replicate[0],b_size=replicate[1],c_size=replicate[2])

        with open(ofname, 'w') as f:
            system.dump('atom_data', f=f, atom_style='atomic')
        natoms = system.natoms
        pafi_info = np.zeros((natoms,10))
        pafi_info[:,0] = np.arange(1, natoms+1) # id: 1 to N
        pafi_info[:,1:4] = system.atoms.pos # X0
        pafi_info[:,4:7] = system.atoms.dX   
        normsq = (pafi_info[:, 4:7]**2).sum() 
        pafi_info[:,4:7] = pafi_info[:,4:7] / np.sqrt(normsq) # dX0/|dX0|
        pafi_info[:,7:10] = system.atoms.ddX / normsq # ddX0/|dX0|^2
        with open(ofname, 'a') as f:
            fmt = '%u'+' %.13g'*9
            np.savetxt(f, pafi_info, fmt=fmt, delimiter=' ',
                       header='\nPafiPath\n', comments='')

    def save(self, fname='path_maker.pkl'):
        '''
        Save PathMaker to file
        '''
        with open(fname, 'wb') as f:
            pickle.dump(self, f)

    def load(self, fname='path_maker.pkl'):
        '''
        Load PathMaker from file
        '''
        with open(fname, 'rb') as f:
            pm = pickle.load(f)
        self.box = pm.box
        self.pbc = pm.pbc
        self.atype = pm.atype
        self.natoms = pm.natoms
        self.path_interpolator = pm.path_interpolator
        self.path_fnames = pm.path_fnames
        self.spline_k = pm.spline_k
        self.reaction_coord = pm.reaction_coord
        self.fitting_quality = pm.fitting_quality


def test_vac():
    import os
    os.chdir('/gauss14/home/cityu/songweili/workspace/research/tial/sma-titanium-aluminide/src/lammps/calc/PAFI/test_vac')

    fnames = './neb_path/dump.neb.*'
    pm = PathMaker()
    pm.fit_path_interpolator(fnames,ftype='dump',check_fit_error=True)
    pm.save('./neb_path/path_maker.pkl')

    pm = PathMaker('./neb_path/path_maker.pkl')
    for r in np.linspace(0, 1, 11):
        pm.make_pafiv1_path(r, odir='./pafiv1_path/')
        pm.make_pafi_path(r, odir='./pafi_path/')

def test_gl():
    import os
    os.chdir('/gauss14/home/cityu/songweili/workspace/research/tial/sma-titanium-aluminide/src/lammps/calc/PAFI/test_gl')

    fnames = './neb_path/dump.neb.*'
    pm = PathMaker()
    pm.fit_path_interpolator(fnames,ftype='dump',check_fit_error=True)
    pm.save('./neb_path/path_maker.pkl')

    pm = PathMaker('./neb_path/path_maker.pkl')
    r_all = np.linspace(0, 1, 11)
    for r in r_all:
        pm.make_pafiv1_path(r, odir='./pafiv1_path/')
        pm.make_pafi_path(r, odir='./pafi_path/')

def main():
    # test_vac()
    # test_gl()

    import argparse, os
    from pathlib import Path
    parser=argparse.ArgumentParser(
        prog="python3 path_maker.py",
        description="Create interpolated path for PAFI calculation",
        epilog="Written by LI Songwei"
    )
    parser.add_argument('--neb_path',type=str,required=True,help='NEB path')
    parser.add_argument('--pickle_file',type=str,required=True,help='Pickle file name')
    parser.add_argument('--reaction_coordinate',type=float,nargs='*',help='Reaction coordinate for output files')
    parser.add_argument('--output_dir',type=str,help='Output PAFI data file directory')
    args=parser.parse_args()
    if os.path.exists(args.pickle_file):
        print("Using existing pickle file:")
        pm=PathMaker(args.pickle_file)
    else:
        pm=PathMaker()
        pm.fit_path_interpolator(fnames=args.neb_path,check_fit_error=True)
        pm.save(args.pickle_file)
    for rc in args.reaction_coordinate:
        rc=float(f'{rc:.3f}')
        pm.make_pafiv1_path(r=rc,odir=args.output_dir)

if __name__ == '__main__':
    main()
