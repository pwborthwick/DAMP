from scipy.linalg import block_diag as kron_diag
import numpy as np

'''
Example of using pySCF integrals with Hugenholz generated code 
'''

#get the pyscf integrals and SCF and MP2 energies
from pyscf import gto, scf, ao2mo

def pyscf_properties():
    mol = gto.Mole(atom  = 'H 0 0 0; F 0 0 1.1',
                  basis = 'ccpvdz',
                  verbose=0)

    mf = mol.RHF().run()
    scf_energy = mf.kernel()
    mp = mf.MP2()
    mp2_energy, _ = mp.kernel()
    
    mo_energy = np.concatenate([mf.mo_energy, mf.mo_energy])
    mo_coeff  = kron_diag(mf.mo_coeff, mf.mo_coeff)
    arg = np.argsort(mo_energy)
    
    eps, mo_coeff = mo_energy[arg], mo_coeff[:, arg]
    
    g = ao2mo.addons.restore(1, mf._eri, mol.nao)
    g = np.kron(np.eye(2), np.kron(np.eye(2), g).transpose())
    g = np.einsum("pqrs,pi,qj,rk,sl->ikjl", g, mo_coeff, mo_coeff,
                                         mo_coeff, mo_coeff, optimize=True)
    g = g - g.transpose(0,1,3,2)

    nocc = mol.nelectron
    
    return eps, g, nocc, scf_energy, mp2_energy

eps, g, nocc, scf_energy, mp2_energy = pyscf_properties()

#hugenholtz module
from outputs.mp2 import damp as mp2

mp2_hugenholtz = mp2(eps, g, nocc)

#results
print('*************************\n* Comparison with pySCF *')
print('*************************')

print('pySCF energies    {:>18.14f} (mp2)   {:>15.10f}'.format(mp2_energy, scf_energy))
print('Hugenholtz energy {:>18.14f} (mp2)'.format(mp2_hugenholtz))

'''
*************************
* Comparison with pySCF *
*************************
pySCF energies     -0.21136746031006 (mp2)    -99.9873974403
Hugenholtz energy  -0.21136746031006 (mp2)
'''
