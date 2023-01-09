def damp(eps, g, nocc):
    #moller-plesset 2

    from numpy import newaxis, reciprocal, einsum

    #orbital slices
    o, v, n = slice(None, nocc), slice(nocc, None), newaxis
    #orbital energy differences
    deltas = {}
    deltas['4'] = reciprocal(eps[o,n,n,n]+eps[n,o,n,n]-eps[n,n,v,n]-eps[n,n,n,v])

    #loop over diagrams
    mp = 0.0
    mp += (pow(2,-2)*einsum('abij,ijab,ijab->',
               g[v,v,o,o], g[o,o,v,v], deltas['4'],
               optimize=True))
    return mp