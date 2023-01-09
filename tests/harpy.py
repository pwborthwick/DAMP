from source.hugenholtz import HUGENHOLTZ
import numpy as np

h2o_sto3g = np.load('h2o-sto3g-hirata.npz')
eps, g, nocc = h2o_sto3g['eps'], h2o_sto3g['g'], h2o_sto3g['nocc']

print('********************************\n* Tests with Hirata Benchmarks *')
print('********************************\n')

print(' theory      correction   %FCI        Hirata')
print('-----------------------------------------------')

fci = -0.0502556615

# mp2
import outputs.mp2 as moller_plesset
mp2 = moller_plesset.damp(eps, g, nocc)
print('   mp2     {:<14.10f}  {:2}%    -0.0358672469'.format(mp2, int(mp2/fci*100) ))

# mp3
import outputs.mp3 as moller_plesset
mp3 = moller_plesset.damp(eps, g, nocc)
print('   mp3     {:<14.10f}  {:2}%    -0.0098015863'.format(mp3, int((mp2+mp3)/fci*100) ))

# mp4
import outputs.mp4 as moller_plesset
mp4 = moller_plesset.damp(eps, g, nocc)
print('   mp4     {:<14.10f}  {:2}%    -0.0030104405'.format(mp4, int((mp2+mp3+mp4)/fci*100) ))

cumulative = mp2 + mp3 + mp4
print('\n     cumulative correction {:>16.10f}'.format(cumulative))
print('     corrected SCF energy    {:<16.10f}'.format(cumulative + -74.9626630861))

'''
********************************
* Tests with Hirata Benchmarks *
********************************

 theory      correction   %FCI        Hirata
-----------------------------------------------
   mp2     -0.0358672477   71%    -0.0358672469
   mp3     -0.0098015864   90%    -0.0098015863
   mp4     -0.0030104405   96%    -0.0030104405

     cumulative correction    -0.0486792746
     corrected SCF energy    -75.0113423607
'''
