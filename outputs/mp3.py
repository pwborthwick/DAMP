def damp(eps, g, nocc):
    #moller-plesset 3

    from numpy import newaxis, reciprocal, einsum

    #orbital slices
    o, v, n = slice(None, nocc), slice(nocc, None), newaxis
    #orbital energy differences
    deltas = {}
    deltas['4'] = reciprocal(eps[o,n,n,n]+eps[n,o,n,n]-eps[n,n,v,n]-eps[n,n,n,v])

    #loop over diagrams
    mp = 0.0
    mp += (pow(2,-3)*einsum('abij,ijkl,klab,ijab,klab->',
               g[v,v,o,o], g[o,o,o,o], g[o,o,v,v], deltas['4'], deltas['4'],
               optimize=True))
    mp += (einsum('abij,icak,jkbc,ijab,jkbc->',
               g[v,v,o,o], g[o,v,v,o], g[o,o,v,v], deltas['4'], deltas['4'],
               optimize=True))
    mp += (pow(2,-3)*einsum('abij,cdab,ijcd,ijab,ijcd->',
               g[v,v,o,o], g[v,v,v,v], g[o,o,v,v], deltas['4'], deltas['4'],
               optimize=True))
    return mp