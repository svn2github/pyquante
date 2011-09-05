"""\
PyQuante2 contains a more structured interface to all of the functions
in PyQuante.

solver = SCF(molecule,**options)

Create a solver that can perform a HF calculation on *molecule*.

General options

Option        Value   Description
--------      -----   -----------
method        HF      Use the HF method for the calculation
              UHF     Use the UHF method for the calculation
              DFT     Use the DFT method for the calculation
              MINDO3  Use the MINDO3 method for the calculation
              UMINDO3 Use the UMINDO3 method for the calculation
              
bfs           None    The basis functions to use. List of CGBF's
basis_data    None    The basis data to use to construct bfs
basis         None    The name of a basis set, e.g. '6-31g**',
                      'sto-3g','cc-pVTZ'
integrals     None    The one- and two-electron integrals to use
                      If not None, S,h,Ints
orbs          None    If not none, the guess orbitals

Options passed into solver.iterate(**options):

Options:      Value   Description
--------      -----   -----------
etol          1e-5    Energy convergence criteria
max_iter      50      Maximum SCF iterations
(do_averaging  True    Use DIIS for accelerated convergence (default)
              False   No convergence acceleration)
etemp         False   Use etemp value for finite temperature DFT (default)
              float   Use (float) for the electron temperature

The test suite at the bottom of the file has examples of usage.
"""

import settings
import unittest,logging

logging.basicConfig(level=settings.LoggingLevel,
                    filename=settings.LoggingFilename)
class SCFIterator:
    def __init__(self,**kwargs):
        self.energy_history = [0]
        self.converged = False
        return

    def iterate(self,ham,**kwargs):
        self.max_iter = kwargs.get('max_iter',settings.MaxIter)
        for self.iter in xrange(1,self.max_iter+1):
            ham.update(**kwargs)
            logging.debug("%d %f" % (self.iter,ham.energy))
            energy_var=abs(ham.energy - self.energy_history[-1]) 
            logging.info("Iteration: %d    Energy: %f    EnergyVar: %f"%
                        (self.iter,ham.energy,energy_var))
            if self.is_converged(ham): break
        if self.iter < self.max_iter:
            logging.info("PyQuante converged in %d iterations" % self.iter)
        else:
            logging.warning("PyQuante failed to converge after %d iterations"
                            % self.max_iter)
        return

    def is_converged(self,ham,**kwargs):
        self.energy = ham.get_energy()
        etol = kwargs.get('etol',settings.ConvergenceCriteria)
        if not self.energy_history:
            self.energy_history.append(self.energy)
            return False
        dE = abs(self.energy-self.energy_history[-1])
        self.converged = dE < etol
        #print self.energy,dE,etol,self.converged
        self.energy_history.append(self.energy)
        return self.converged

    def __repr__(self):
        lstr = ["Iterator information:"]
        lstr.extend([str(en) for en in self.energy_history])
        if self.converged:
            lstr.append("The iterator is converged")
        else:
            lstr.append("The iterator is not converged")
        return "\n".join(lstr)

class Integrals:
    def __init__(self,molecule,basis_set,**kwargs):
        from PyQuante.Ints import getints
        integrals = kwargs.get("integrals")
        nbf = len(basis_set)
        self.S, self.h, self.ERI = getints(basis_set.get(),molecule,**kwargs)
        return
    

    def get(self): return self.S, self.h, self.ERI
    def get_S(self): return self.S
    def get_h(self): return self.h
    def get_ERI(self): return self.ERI

class BasisSet:
    def __init__(self,molecule,**kwargs):
        from PyQuante.Ints import getbasis
        from PyQuante.Basis.Tools import get_basis_data

        self.bfs = getbasis(molecule,**kwargs)
        logging.info("%d basis functions" % len(self.bfs))
        return
    def __repr__(self): return 'Gaussian basis set with %d bfns' %  len(self.bfs)
    def __len__(self): return len(self.bfs)

    def get(self): return self.bfs


########## Hamiltonian ##########

def HamiltonianFactory(molecule,**kwargs):
    method = kwargs.get('method',settings.HamiltonianMethod)
    if method == "UHF":
        return UHFHamiltonian(molecule,**kwargs)
    elif method == 'ROHF':
        return ROHFHamiltonian(molecule,**kwargs)
    elif method == "DFT":
        return DFTHamiltonian(molecule,**kwargs)
    elif method == 'MINDO3':
        return MINDO3Hamiltonian(molecule,**kwargs)
    elif method == 'UMINDO3':
        return UMINDO3Hamiltonian(molecule,**kwargs)
    return HFHamiltonian(molecule,**kwargs)
# Convenience function
def SCF(molecule,**kwargs): return HamiltonianFactory(molecule,**kwargs)

class AbstractHamiltonian:
    def __init__(self,molecule,**kwargs):
        raise Exception("AbstractHamiltonian::__init__")
    def update(self,**kwargs):
        raise Exception("AbstractHamiltonian::update")
    def iterate(self,**kwargs):
        raise Exception("AbstractHamiltonian::iterate")
    def get_energy(self,**kwargs):
        raise Exception("AbstractHamiltonian::get_energy")

class HFHamiltonian(AbstractHamiltonian):
    method='HF'
    def __init__(self,molecule,**kwargs):
        from PyQuante.Convergence import DIIS
        self.molecule = molecule
        logging.info("HF calculation on system %s" % self.molecule.name)
        self.basis_set = BasisSet(molecule,**kwargs)
        self.integrals = Integrals(molecule,self.basis_set,**kwargs)
        self.iterator = SCFIterator()
        self.h = self.integrals.get_h()
        self.S = self.integrals.get_S()
        self.ERI = self.integrals.get_ERI()
        self.Enuke = molecule.get_enuke()
        self.F = self.h
        self.dmat = None
        self.entropy = None
        self.DoAveraging = kwargs.get('DoAveraging',settings.Averaging)
        if self.DoAveraging:
            self.Averager = DIIS(self.S)
        nel = molecule.get_nel()
        nclosed,nopen = molecule.get_closedopen()
        logging.info("Nclosed/open = %d, %d" % (nclosed,nopen))
        self.solver = SolverFactory(nel,nclosed,nopen,self.S,**kwargs)
        return

    def __repr__(self):
        lstr = ['Hamiltonian constructed for method %s' % self.method,
                repr(self.molecule),
                repr(self.basis_set),
                repr(self.iterator)]
        return '\n'.join(lstr)

    def get_energy(self): return self.energy
    def iterate(self,**kwargs):
        self.iterator.iterate(self,**kwargs)
        if self.iterator.converged:
            logging.info("Final HF energy for system %s is %f"%(self.molecule.name,self.energy))
        return

    def update(self,**kwargs):
        from PyQuante.LA2 import trace2
        from PyQuante.Ints import getJ,getK

        if self.DoAveraging and self.dmat is not None:
            self.F = self.Averager.getF(self.F,self.dmat)
        self.dmat,self.entropy = self.solver.solve(self.F,**kwargs)
        D = self.dmat
        
        self.J = getJ(self.ERI,D)
        self.Ej = 2*trace2(D,self.J)
        self.K = getK(self.ERI,D)
        self.Exc = -trace2(D,self.K)
        self.Eone = 2*trace2(D,self.h)
        self.F = self.h + 2*self.J - self.K
        self.energy = self.Eone + self.Ej + self.Exc + self.Enuke + self.entropy
        return

class DFTHamiltonian(AbstractHamiltonian):
    method='DFT'
    def __init__(self,molecule,**kwargs):
        from PyQuante.Convergence import DIIS
        from PyQuante.DFunctionals import need_gradients
        self.molecule = molecule
        logging.info("DFT calculation on system %s" % self.molecule.name)
        self.basis_set = BasisSet(molecule,**kwargs)
        self.integrals = Integrals(molecule,self.basis_set,**kwargs)
        self.iterator = SCFIterator()
        self.h = self.integrals.get_h()
        self.S = self.integrals.get_S()
        self.ERI = self.integrals.get_ERI()
        self.Enuke = molecule.get_enuke()
        self.nel = molecule.get_nel()
        self.F = self.h
        self.functional = kwargs.get('functional',settings.DFTFunctional)
        kwargs['do_grad_dens'] = need_gradients[self.functional]
        self.setup_grid(molecule,self.basis_set.get(),**kwargs)
        self.dmat = None
        self.entropy = None
        self.DoAveraging = kwargs.get('DoAveraging',settings.DFTAveraging)
        if self.DoAveraging:
            self.Averager = DIIS(self.S)
        nel = molecule.get_nel()
        nclosed,nopen = molecule.get_closedopen()
        logging.info("Nclosed/open = %d, %d" % (nclosed,nopen))
        self.solver = SolverFactory(nel,nclosed,nopen,self.S,**kwargs)
        return
        
    def __repr__(self):
        lstr = ['Hamiltonian constructed for method %s' % self.method,
                repr(self.molecule),
                repr(self.basis_set),
                repr(self.iterator)]
        return '\n'.join(lstr)

    def get_energy(self): return self.energy
    def iterate(self,**kwargs): return self.iterator.iterate(self,**kwargs)

    def setup_grid(self,molecule,bfs,**kwargs):
        #from PyQuante.MolecularGrid import MolecularGrid
        from PyQuante.MG2 import MG2 as MolecularGrid
        grid_nrad = kwargs.get('grid_nrad',settings.DFTGridRadii)
        grid_fineness = kwargs.get('grid_fineness',settings.DFTGridFineness)
        self.gr = MolecularGrid(molecule,grid_nrad,grid_fineness,**kwargs) 
        self.gr.set_bf_amps(bfs)
        return

    def update(self,**kwargs):
        from PyQuante.LA2 import trace2
        from PyQuante.Ints import getJ
        from PyQuante.dft import getXC

        #self.DoAveraging = kwargs.get('DoAveraging',True)
        #if self.DoAveraging:
        #    self.Averager = DIIS(self.S)

        if self.DoAveraging and self.dmat is not None:
            self.F = self.Averager.getF(self.F,self.dmat)
        self.dmat,self.entropy = self.solver.solve(self.F,**kwargs)
        D = self.dmat
        
        self.gr.setdens(D)
        self.J = getJ(self.ERI,D)
        self.Ej = 2*trace2(D,self.J)

        self.Exc,self.XC = getXC(self.gr,self.nel,
                                 functional=self.functional)

        self.Eone = 2*trace2(D,self.h)

        self.F = self.h+2*self.J+self.XC
        self.energy = self.Eone + self.Ej + self.Exc + self.Enuke + self.entropy
        return

class UHFHamiltonian(AbstractHamiltonian):
    method='UHF'
    def __init__(self,molecule,**kwargs):
        self.molecule = molecule
        logging.info("UHF calculation on system %s" % self.molecule.name)
        self.basis_set = BasisSet(molecule,**kwargs)
        self.integrals = Integrals(molecule,self.basis_set,**kwargs)
        self.iterator = SCFIterator()
        self.h = self.integrals.get_h()
        self.S = self.integrals.get_S()
        self.ERI = self.integrals.get_ERI()
        self.Enuke = molecule.get_enuke()

        self.Fa = self.h
        self.Fb = self.h

        self.amat = None
        self.bmat = None
        self.entropy = None
        nalpha,nbeta = molecule.get_alphabeta()
        logging.info("Nalpha/beta = %d, %d" % (nalpha,nbeta))
        self.solvera = SolverFactory(2*nalpha,nalpha,0,self.S,**kwargs)
        self.solverb = SolverFactory(2*nbeta,nbeta,0,self.S,**kwargs)
        return

    def __repr__(self):
        lstr = ['Hamiltonian constructed for method %s' % self.method,
                repr(self.molecule),
                repr(self.basis_set),
                repr(self.iterator)]
        return '\n'.join(lstr)

    def get_energy(self): return self.energy
    def iterate(self,**kwargs): return self.iterator.iterate(self,**kwargs)

    def update(self,**kwargs):
        from PyQuante.LA2 import trace2
        from PyQuante.Ints import getJ,getK

        self.amat,entropya = self.solvera.solve(self.Fa)
        self.bmat,entropyb = self.solverb.solve(self.Fb)

        Da = self.amat
        Db = self.bmat

        D = Da+Db
        self.entropy = 0.5*(entropya+entropyb)

        self.J = getJ(self.ERI,D)
        self.Ej = 0.5*trace2(D,self.J)
        self.Ka = getK(self.ERI,Da)
        self.Kb = getK(self.ERI,Db)
        self.Exc = -0.5*(trace2(Da,self.Ka)+trace2(Db,self.Kb))
        self.Eone = trace2(D,self.h)
        self.Fa = self.h + self.J - self.Ka
        self.Fb = self.h + self.J - self.Kb
        self.energy = self.Eone + self.Ej + self.Exc + self.Enuke + self.entropy
        return

class ROHFHamiltonian(AbstractHamiltonian):
    method='ROHF'
    def __init__(self,molecule,**kwargs):
        self.molecule = molecule
        logging.info("ROHF calculation on system %s" % self.molecule.name)
        self.basis_set = BasisSet(molecule,**kwargs)
        self.integrals = Integrals(molecule,self.basis_set,**kwargs)
        self.iterator = SCFIterator()
        self.h = self.integrals.get_h()
        self.S = self.integrals.get_S()
        self.ERI = self.integrals.get_ERI()
        self.Enuke = molecule.get_enuke()

        self.orbs = None
        self.norbs = len(self.basis_set)

        self.nalpha,self.nbeta = molecule.get_alphabeta()
        logging.info("Nalpha/beta = %d, %d" % (self.nalpha,self.nbeta))
        return

    def __repr__(self):
        lstr = ['Hamiltonian constructed for method %s' % self.method,
                repr(self.molecule),
                repr(self.basis_set),
                repr(self.iterator)]
        return '\n'.join(lstr)

    def get_energy(self): return self.energy
    def iterate(self,**kwargs): return self.iterator.iterate(self,**kwargs)

    def update(self,**kwargs):
        from PyQuante.Ints import getJ,getK
        from PyQuante.LA2 import geigh,mkdens
        from PyQuante.rohf import ao2mo
        from PyQuante.hartree_fock import get_energy
        from PyQuante.NumWrap import eigh,matrixmultiply

        if self.orbs is None:
            self.orbe,self.orbs = geigh(self.h, self.S)
        Da = mkdens(self.orbs,0,self.nalpha)
        Db = mkdens(self.orbs,0,self.nbeta)

        Ja = getJ(self.ERI,Da)
        Jb = getJ(self.ERI,Db)
        Ka = getK(self.ERI,Da)
        Kb = getK(self.ERI,Db)
        Fa = self.h+Ja+Jb-Ka
        Fb = self.h+Ja+Jb-Kb
        energya = get_energy(self.h,Fa,Da)
        energyb = get_energy(self.h,Fb,Db)
        self.energy = (energya+energyb)/2 + self.Enuke

        Fa = ao2mo(Fa,self.orbs)
        Fb = ao2mo(Fb,self.orbs)

        # Building the approximate Fock matrices in the MO basis
        F = 0.5*(Fa+Fb)
        K = Fb-Fa

        # The Fock matrix now looks like
        #      F-K    |  F + K/2  |    F
        #   ---------------------------------
        #    F + K/2  |     F     |  F - K/2
        #   ---------------------------------
        #       F     |  F - K/2  |  F + K

        # Make explicit slice objects to simplify this
        do = slice(0,self.nbeta)
        so = slice(self.nbeta,self.nalpha)
        uo = slice(self.nalpha,self.norbs)
        F[do,do] -= K[do,do]
        F[uo,uo] += K[uo,uo]
        F[do,so] += 0.5*K[do,so]
        F[so,do] += 0.5*K[so,do]
        F[so,uo] -= 0.5*K[so,uo]
        F[uo,so] -= 0.5*K[uo,so]

        self.orbe,mo_orbs = eigh(F)
        self.orbs = matrixmultiply(self.orbs,mo_orbs)
        
        return

class MINDO3Hamiltonian(AbstractHamiltonian):
    method='MINDO3'
    def __init__(self,molecule,**kwargs):
        from PyQuante.MINDO3 import initialize, get_nbf, get_reference_energy,\
             get_F0, get_nel,get_open_closed,get_enuke,get_guess_D
        self.molecule = molecule
        logging.info("MINDO3 calculation on system %s" % self.molecule.name)
        self.iterator = SCFIterator()
        self.charge = self.molecule.charge
        self.multiplicity = self.molecule.multiplicity
        # This is an ugly-ish hack to deal with the brain-dead
        #  way that MINDO3.get_open_closed works
        if self.multiplicity == 1: self.multiplicity=None
        # Ultimately I should subclass Atom for MINDO3Atom
        self.molecule = initialize(self.molecule)
        self.nel = get_nel(self.molecule,self.charge)
        self.nclosed,self.nopen = get_open_closed(self.nel,self.multiplicity)
        self.Enuke = get_enuke(self.molecule)
        self.energy = 0
        self.method = "MINDO3"
        self.nbf = get_nbf(self.molecule)
        self.eref = get_reference_energy(self.molecule)
        self.F0 = get_F0(self.molecule)
        self.F = self.F0
        self.D = get_guess_D(self.molecule)
        logging.info("Nel = %d Nclosed = %d Nopen = %d Enuke = %f Nbf = %d"
                     % (self.nel,self.nclosed,self.nopen,self.Enuke,self.nbf))
        return

    def __repr__(self):
        lstr = ['Hamiltonian constructed for method %s' % self.method,
                repr(self.molecule),
                'Implicit MINDO3 basis set with %d bfns' % self.nbf,
                repr(self.iterator)]
        return '\n'.join(lstr)

    def get_energy(self): return self.energy
    def iterate(self,**kwargs): return self.iterator.iterate(self,**kwargs)

    def update(self,**kwargs):
        self.update_fock()
        self.calculate_energy()
        self.solve_fock()
        self.update_density()

    def update_fock(self):
        from PyQuante.MINDO3 import get_F1, get_F2
        avg = 0.25
        Fold = self.F
        self.F1 = get_F1(self.molecule,self.D)
        self.F2 = get_F2(self.molecule,self.D)
        self.F = self.F0+self.F1+self.F2
        #self.F = avg*self.F + (1-avg)*Fold
        
    def solve_fock(self):
        from PyQuante.NumWrap import eigh
        self.orbe,self.orbs = eigh(self.F)

    def update_density(self):
        from PyQuante.LA2 import mkdens
        self.D = 2*mkdens(self.orbs,0,self.nclosed)
        
    def calculate_energy(self):
        from PyQuante.LA2 import trace2
        from PyQuante.MINDO3 import ev2kcal
        self.Eel = 0.5*trace2(self.D,self.F0+self.F)
        self.Etot = self.Eel+self.Enuke
        self.energy = self.Etot*ev2kcal+self.eref
        
class UMINDO3Hamiltonian(AbstractHamiltonian):
    def __init__(self,molecule,**kwargs):
        from PyQuante.MINDO3 import initialize, get_nbf, get_reference_energy,\
             get_F0, get_nel,get_open_closed,get_enuke,get_guess_D
        self.molecule = molecule
        logging.info("uMINDO3 calculation on system %s" % self.molecule.name)
        self.iterator = SCFIterator()
        self.charge = self.molecule.charge
        self.multiplicity = self.molecule.multiplicity
        # This is an ugly-ish hack to deal with the brain-dead
        #  way that MINDO3.get_open_closed works
        if self.multiplicity == 1: self.multiplicity=None
        # Ultimately I should subclass Atom for MINDO3Atom
        self.molecule = initialize(self.molecule)
        self.nel = get_nel(self.molecule,self.charge)
        self.nclosed,self.nopen = get_open_closed(self.nel,self.multiplicity)
        logging.info("Nclosed/open = %d, %d" % (self.nclosed,self.nopen))
        self.Enuke = get_enuke(self.molecule)
        self.energy = 0
        self.method = "MINDO3"
        self.nbf = get_nbf(self.molecule)
        self.eref = get_reference_energy(self.molecule)
        self.F0 = get_F0(self.molecule)
        self.Fa = self.Fb = self.F0
        self.nalpha = self.nclosed+self.nopen
        self.nbeta = self.nclosed
        self.Da = self.Db = 0.5*get_guess_D(self.molecule)
        self.start = True
        return

    def __repr__(self):
        lstr = ['Hamiltonian constructed for method %s' % self.method,
                repr(self.molecule),
                'Implicit MINDO3 basis set with %d bfns' % self.nbf,
                repr(self.iterator)]
        return '\n'.join(lstr)

    def get_energy(self): return self.energy
    def iterate(self,**kwargs): return self.iterator.iterate(self,**kwargs)

    def update(self,**kwargs):
        self.solve_fock()
        self.update_density()
        self.update_fock()
        self.calculate_energy()

    def update_fock(self):
        from PyQuante.MINDO3 import get_F1_open, get_F2_open
        F1a = get_F1_open(self.molecule,self.Da,self.Db)
        F1b = get_F1_open(self.molecule,self.Db,self.Da)
        F2a = get_F2_open(self.molecule,self.Da,self.Db)
        F2b = get_F2_open(self.molecule,self.Db,self.Da)
        self.Fa = self.F0+F1a+F2a
        self.Fb = self.F0+F1b+F2b
        return

    def solve_fock(self):
        from PyQuante.NumWrap import eigh
        from PyQuante.LA2 import mkdens
        self.orbea,self.orbsa = eigh(self.Fa)
        self.orbeb,self.orbsb = eigh(self.Fb)

    def update_density(self):
        from PyQuante.LA2 import mkdens
        if self.start:
            self.start = False
        else:
            self.Da = mkdens(self.orbsa,0,self.nalpha)
            self.Db = mkdens(self.orbsb,0,self.nbeta)

    def calculate_energy(self):
        from PyQuante.LA2 import trace2
        from PyQuante.MINDO3 import ev2kcal
        self.Eel = 0.5*trace2(self.Da,self.F0+self.Fa)+\
                   0.5*trace2(self.Db,self.F0+self.Fb)
        self.Etot = self.Eel+self.Enuke
        self.energy = self.Etot*ev2kcal+self.eref

########## Solver ##########

def SolverFactory(nel,nclosed,nopen,S,**kwargs):
    if kwargs.get("SolverConstructor"):
        # We can override all of this and pass in an explicit solver constructor:
        return kwargs["SolverConstructor"](nel,nclosed,nopen,S,**kwargs)
    if kwargs.get('etemp',settings.ElectronTemperature):
        return FermiDiracSolver(nel,nclosed,nopen,S,**kwargs)
    return BasicSolver(nel,nclosed,nopen,S,**kwargs)

class AbstractSolver:
    def __init__(self,S,**kwargs):
        raise Exception("AbstractSolver::__init__")
    def solve(self,ham,**kwargs):
        raise Exception("AbstractSolver::solve")

class BasicSolver(AbstractSolver):
    def __init__(self,nel,nclosed,nopen,S,**kwargs):
        self.S = S
        self.nel = nel
        self.nclosed = nclosed
        self.nopen = nopen
        return

    def solve(self,H,**kwargs):
        from PyQuante.LA2 import geigh,mkdens_spinavg
        self.orbe,self.orbs = geigh(H,self.S)
        self.D = mkdens_spinavg(self.orbs,self.nclosed,self.nopen)
        self.entropy = 0
        return self.D,self.entropy

class FermiDiracSolver(AbstractSolver):
    def __init__(self,nel,nclosed,nopen,S,**kwargs):
        self.S = S
        self.nel = nel
        self.nclosed = nclosed
        self.nopen = nopen
        self.etemp = kwargs.get('etemp',settings.ElectronTemperature)
        return

    def solve(self,H,**kwargs):
        from PyQuante.LA2 import geigh
        from PyQuante.fermi_dirac import mkdens_fermi
        self.orbe,self.orbs = geigh(H,self.S)
        self.D,self.entropy = mkdens_fermi(self.nel,self.orbe,
                                           self.orbs,self.etemp)
        return self.D,self.entropy

class SubspaceSolver(AbstractSolver):
    def __init__(self,nel,nclosed,nopen,S,**kwargs):
        self.S = S
        self.nel = nel
        self.nclosed = nclosed
        self.nopen = nopen
        self.first_iteration = True
        self.solver = kwargs.get("solver")

        # Determine nroots, which Davidson (and perhaps others) needs:
        self.pass_nroots = kwargs.get("pass_nroots")
        self.nvirt = kwargs.get("nvirt",settings.SubspaceVirtualOrbs) # solve for 1 virtual orbital by default
        self.nroots = self.nclosed + self.nopen + self.nvirt

        if not self.solver:
            from PyQuante.NumWrap import eigh
            self.solver = eigh
        return

    def solve(self,H,**kwargs):
        from PyQuante.LA2 import mkdens_spinavg,simx,geigh
        from PyQuante.NumWrap import matrixmultiply,eigh
        if self.first_iteration:
            self.first_iteration = False
            self.orbe,self.orbs = geigh(H,self.S)
        else:
            Ht = simx(H,self.orbs)
            if self.pass_nroots:
                self.orbe,orbs = self.solver(Ht,self.nroots)
            else:
                self.orbe,orbs = self.solver(Ht)
            self.orbs = matrixmultiply(self.orbs,orbs)
        self.D = mkdens_spinavg(self.orbs,self.nclosed,self.nopen)
        self.entropy = 0
        return self.D,self.entropy

class DmatSolver(AbstractSolver):
    def __init__(self,nel,nclosed,nopen,S,**kwargs):
        self.S = S
        self.nel = nel
        self.nclosed = nclosed
        self.nopen = nopen
        self.solver = kwargs.get("solver")
        if not self.solver:
            from PyQuante.DMP import TCP
            self.solver = TCP
        return

    def solve(self,H,**kwargs):
        solver = self.solver(H,self.nclosed,self.S)
        solver.iterate()
        self.entropy = 0
        self.D = solver.D
        return self.D,self.entropy

class UnitTests(unittest.TestCase):
    def setUp(self):
        from PyQuante.Molecule import Molecule
        self.h2 = Molecule('H2',atomlist=[(1,(0.35,0,0)),(1,(-0.35,0,0))],
                           units='Angs')
        self.he = Molecule('He',atomlist = [(2,(0,0,0))])
        self.li = Molecule('Li',atomlist = [(3,(0,0,0))],multiplicity=2)
        self.li_p = Molecule('Li+',atomlist = [(3,(0,0,0))],charge=1)
        self.li_m = Molecule('Li-',atomlist = [(3,(0,0,0))],charge=-1)
        self.h2o = Molecule('h2o',[(8,(0,0,0)),(1,(1.,0,0)),(1,(0,1.,0))],
                            units="Angstrom")
        self.oh = Molecule('oh',[(8,(0,0,0)),(1,(1.,0,0))],
                           units="Angstrom")

    def testH2HF(self):
        h2_hf = SCF(self.h2,method='HF')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def testH2HFFT(self):
        h2_hf = SCF(self.h2,method='HF',etemp=1e4)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130502,4)

    def testHeHF(self):
        he_hf = SCF(self.he,method='HF')
        he_hf.iterate()
        self.assertAlmostEqual(he_hf.energy,-2.855260,3)

    def testLiHF(self):
        li_hf = SCF(self.li,method='HF')
        li_hf.iterate()
        return

    def testLipHF(self):
        li_p_hf = SCF(self.li_p,method='HF')
        li_p_hf.iterate()
        self.assertAlmostEqual(li_p_hf.energy,-7.235536,4)

    def testLimHF(self):
        li_m_hf = SCF(self.li_m,method='HF')
        li_m_hf.iterate()
        self.assertAlmostEqual(li_m_hf.energy,-7.407030,4)

    def testH2LDA(self):
        h2_lda = SCF(self.h2,method='DFT',functional="SVWN")
        h2_lda.iterate()
        self.assertAlmostEqual(h2_lda.energy,-1.132799,4)

    def testH2LDAFT(self):
        h2_lda = SCF(self.h2,method='DFT',functional="SVWN",etemp=1e4)
        h2_lda.iterate()
        self.assertAlmostEqual(h2_lda.energy,-1.132558,4)

    def testLiLDA(self):
        li_lda = SCF(self.li,method='DFT',functional="SVWN")
        li_lda.iterate()
        self.assertAlmostEqual(li_lda.energy,-7.332050,4)

    def testLiLDAFT(self):
        li_lda = SCF(self.li,method='DFT',functional="SVWN",etemp=1e4)
        li_lda.iterate()
        self.assertAlmostEqual(li_lda.energy,-7.349422,4)

    def testLiUHF(self):
        li_uhf = SCF(self.li,method='UHF')
        li_uhf.iterate()
        self.assertAlmostEqual(li_uhf.energy,-7.431364,4)

    def testLiUHFFT(self):
        li_uhf = SCF(self.li,method="UHF",etemp=1e4)
        li_uhf.iterate()
        # No test, since I don't really know what the energy should be:
        # finite temperature HF is kind of a hack. But this at least
        # tests that the program runs
        return

    ########## Solver tests ##########

    def testSubspaceSolver(self):
        h2_hf = SCF(self.h2,method='HF',SolverConstructor=SubspaceSolver)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)
    
    def testDavidsonSolver(self):
        from PyQuante.Solvers import davidson
        h2_hf = SCF(self.h2,method='HF',SolverConstructor=SubspaceSolver,
                    solver=davidson,pass_nroots=True)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def testJacobiSolver(self):
        from PyQuante.Solvers import jacobi
        h2_hf = SCF(self.h2,method='HF',SolverConstructor=SubspaceSolver,
                    solver=jacobi)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)
    
    def testTCPSolver(self):
        h2_hf = SCF(self.h2,method='HF',SolverConstructor=DmatSolver)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def testTRPSolver(self):
        from PyQuante.DMP import TRP
        h2_hf = SCF(self.h2,method='HF',SolverConstructor=DmatSolver,
                    solver=TRP)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def testCPSolver(self):
        from PyQuante.DMP import CP
        h2_hf = SCF(self.h2,method='HF',SolverConstructor=DmatSolver,
                    solver=CP)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)
    
    def testMCWSolver(self):
        from PyQuante.DMP import McWeeny
        h2_hf = SCF(self.h2,method='HF',SolverConstructor=DmatSolver,
                    solver=McWeeny)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    ########## Basis set tests ##########

    def testSTO3G(self):
        h2_hf = SCF(self.h2,method='HF',basis='sto-3g')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.117349,4)

    def testSTO6G(self):
        h2_hf = SCF(self.h2, method="HF",basis='sto-3g')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.117349,4)

    def test321G(self):
        h2_hf = SCF(self.h2, method="HF",basis='3-21g')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.121998,4)

    def test631Gss(self):
        h2_hf = SCF(self.h2, method="HF",basis='6-31g**')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def test631Gppss(self):
        h2_hf = SCF(self.h2, method="HF",basis='6-31g++**')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130553,4)

    def test631Gdp(self):
        h2_hf = SCF(self.h2, method="HF",basis='6-31G(d,p)')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def testVDZ(self):
        h2_hf = SCF(self.h2, method="HF",basis='cc-pvdz')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.126923,4)

    def testVTZ(self):
        h2_hf = SCF(self.h2, method="HF",basis='cc-pvtz')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.132136,4)

    def testDZVP(self):
        h2_hf = SCF(self.h2, method="HF",basis='dzvp')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.126728,4)

    def test6311G(self):
        h2_hf = SCF(self.h2, method="HF",basis='6-311G**')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.131516,4)

    def test6311Gdp(self):
        h2_hf = SCF(self.h2, method="HF",basis='6-311G++(2d,2p)')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.132122,4)

    def test6311G3d3p(self):
        h2_hf = SCF(self.h2, method="HF",basis='6-311G++(3d,3p)')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.132166,4)

    ########## Misc Tests ##########

    def testH2OMINDO(self):
        h2o_mindo3 = SCF(self.h2o,method="MINDO3")
        h2o_mindo3.iterate()
        self.assertAlmostEqual(h2o_mindo3.energy,-48.826208,2)

    def testOHMINDO(self):
        oh_mindo = SCF(self.oh,method="UMINDO3")
        oh_mindo.iterate()
        self.assertAlmostEqual(oh_mindo.energy,18.1258,2)

def test():
    logging.basicConfig(level=logging.DEBUG,format="%(message)s")
    suite = unittest.TestLoader().loadTestsFromTestCase(UnitTests)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()

if __name__ == '__main__': test()
