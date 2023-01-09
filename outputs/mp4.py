def damp(eps, g, nocc):
    #moller-plesset 4

    from numpy import newaxis, reciprocal, einsum

    #orbital slices
    o, v, n = slice(None, nocc), slice(nocc, None), newaxis
    #orbital energy differences
    deltas = {}
    deltas['2'] = reciprocal(eps[o,n]-eps[n,v])

    deltas['4'] = reciprocal(eps[o,n,n,n]+eps[n,o,n,n]-eps[n,n,v,n]-eps[n,n,n,v])

    deltas['6'] = reciprocal(eps[o,n,n,n,n,n]+eps[n,o,n,n,n,n]+eps[n,n,o,n,n,n]-eps[n,n,n,v,n,n]-eps[n,n,n,n,v,n]-eps[n,n,n,n,n,v])

    deltas['8'] = reciprocal(eps[o,n,n,n,n,n,n,n]+eps[n,o,n,n,n,n,n,n]+eps[n,n,o,n,n,n,n,n]+eps[n,n,n,o,n,n,n,n]-eps[n,n,n,n,v,n,n,n]-eps[n,n,n,n,n,v,n,n]-eps[n,n,n,n,n,n,v,n]-eps[n,n,n,n,n,n,n,v])

    #loop over diagrams
    mp = 0.0
    mp += (-pow(2,-2)*einsum('abij,cdkl,ikcd,jlab,ijab,ijklcdab,jlab->',
               g[v,v,o,o], g[v,v,o,o], g[o,o,v,v], g[o,o,v,v], deltas['4'], deltas['8'], deltas['4'],
               optimize=True))
    mp += (-pow(2,-2)*einsum('abij,cdkl,klac,ijbd,ijab,ijklacbd,ijbd->',
               g[v,v,o,o], g[v,v,o,o], g[o,o,v,v], g[o,o,v,v], deltas['4'], deltas['8'], deltas['4'],
               optimize=True))
    mp += (pow(2,-4)*einsum('abij,cdkl,ijcd,klab,ijab,ijklcdab,klab->',
               g[v,v,o,o], g[v,v,o,o], g[o,o,v,v], g[o,o,v,v], deltas['4'], deltas['8'], deltas['4'],
               optimize=True))
    mp += (einsum('abij,cdkl,ikac,jlbd,ijab,ijklacbd,jlbd->',
               g[v,v,o,o], g[v,v,o,o], g[o,o,v,v], g[o,o,v,v], deltas['4'], deltas['8'], deltas['4'],
               optimize=True))
    mp += (pow(2,-4)*einsum('abij,cdkl,klab,ijcd,ijab,ijklabcd,ijcd->',
               g[v,v,o,o], g[v,v,o,o], g[o,o,v,v], g[o,o,v,v], deltas['4'], deltas['8'], deltas['4'],
               optimize=True))
    mp += (-pow(2,-2)*einsum('abij,cdkl,ijac,klbd,ijab,ijklacbd,klbd->',
               g[v,v,o,o], g[v,v,o,o], g[o,o,v,v], g[o,o,v,v], deltas['4'], deltas['8'], deltas['4'],
               optimize=True))
    mp += (-pow(2,-2)*einsum('abij,cdkl,ikab,jlcd,ijab,ijklabcd,jlcd->',
               g[v,v,o,o], g[v,v,o,o], g[o,o,v,v], g[o,o,v,v], deltas['4'], deltas['8'], deltas['4'],
               optimize=True))
    mp += (pow(2,-2)*einsum('abij,ickl,klcm,jmab,ijab,jklcab,jmab->',
               g[v,v,o,o], g[o,v,o,o], g[o,o,v,o], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (pow(2,-2)*einsum('abij,cdak,kecd,ijbe,ijab,ijkcdb,ijbe->',
               g[v,v,o,o], g[v,v,v,o], g[o,v,v,v], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (pow(2,-4)*einsum('abij,ijkl,klmn,mnab,ijab,klab,mnab->',
               g[v,v,o,o], g[o,o,o,o], g[o,o,o,o], g[o,o,v,v], deltas['4'], deltas['4'], deltas['4'],
               optimize=True))
    mp += (einsum('abij,icak,kdcl,jlbd,ijab,jkcb,jlbd->',
               g[v,v,o,o], g[o,v,v,o], g[o,v,v,o], g[o,o,v,v], deltas['4'], deltas['4'], deltas['4'],
               optimize=True))
    mp += (pow(2,-4)*einsum('abij,cdab,efcd,ijef,ijab,ijcd,ijef->',
               g[v,v,o,o], g[v,v,v,v], g[v,v,v,v], g[o,o,v,v], deltas['4'], deltas['4'], deltas['4'],
               optimize=True))
    mp += (pow(2,-2)*einsum('abij,ijak,kclm,lmbc,ijab,kb,lmbc->',
               g[v,v,o,o], g[o,o,v,o], g[o,v,o,o], g[o,o,v,v], deltas['4'], deltas['2'], deltas['4'],
               optimize=True))
    mp += (pow(2,-2)*einsum('abij,icab,deck,jkde,ijab,jc,jkde->',
               g[v,v,o,o], g[o,v,v,v], g[v,v,v,o], g[o,o,v,v], deltas['4'], deltas['2'], deltas['4'],
               optimize=True))
    mp += (-pow(2,-2)*einsum('abij,ickl,jdab,klcd,ijab,jklabc,klcd->',
               g[v,v,o,o], g[o,v,o,o], g[o,v,v,v], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (-pow(2,-2)*einsum('abij,cdak,ijbl,klcd,ijab,ijkbcd,klcd->',
               g[v,v,o,o], g[v,v,v,o], g[o,o,v,o], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (pow(2,-4)*einsum('abij,ijkl,cdab,klcd,ijab,klab,klcd->',
               g[v,v,o,o], g[o,o,o,o], g[v,v,v,v], g[o,o,v,v], deltas['4'], deltas['4'], deltas['4'],
               optimize=True))
    mp += (einsum('abij,icak,jdbl,klcd,ijab,jkbc,klcd->',
               g[v,v,o,o], g[o,v,v,o], g[o,v,v,o], g[o,o,v,v], deltas['4'], deltas['4'], deltas['4'],
               optimize=True))
    mp += (pow(2,-4)*einsum('abij,cdab,ijkl,klcd,ijab,ijcd,klcd->',
               g[v,v,o,o], g[v,v,v,v], g[o,o,o,o], g[o,o,v,v], deltas['4'], deltas['4'], deltas['4'],
               optimize=True))
    mp += (-pow(2,-2)*einsum('abij,ijak,cdbl,klcd,ijab,kb,klcd->',
               g[v,v,o,o], g[o,o,v,o], g[v,v,v,o], g[o,o,v,v], deltas['4'], deltas['2'], deltas['4'],
               optimize=True))
    mp += (-pow(2,-2)*einsum('abij,icab,jdkl,klcd,ijab,jc,klcd->',
               g[v,v,o,o], g[o,v,v,v], g[o,v,o,o], g[o,o,v,v], deltas['4'], deltas['2'], deltas['4'],
               optimize=True))
    mp += (pow(2,-1)*einsum('abij,ickl,jkcm,lmab,ijab,jklcab,lmab->',
               g[v,v,o,o], g[o,v,o,o], g[o,o,v,o], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (pow(2,-1)*einsum('abij,ickl,klam,jmbc,ijab,jklabc,jmbc->',
               g[v,v,o,o], g[o,v,o,o], g[o,o,v,o], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (-einsum('abij,ickl,kdac,jlbd,ijab,jklacb,jlbd->',
               g[v,v,o,o], g[o,v,o,o], g[o,v,v,v], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (-einsum('abij,cdak,ikcl,jlbd,ijab,ijkcbd,jlbd->',
               g[v,v,o,o], g[v,v,v,o], g[o,o,v,o], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (pow(2,-1)*einsum('abij,cdak,iecd,jkbe,ijab,ijkcdb,jkbe->',
               g[v,v,o,o], g[v,v,v,o], g[o,v,v,v], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (pow(2,-1)*einsum('abij,cdak,kebc,ijde,ijab,ijkbcd,ijde->',
               g[v,v,o,o], g[v,v,v,o], g[o,v,v,v], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (einsum('abij,ickl,jkam,lmbc,ijab,jklabc,lmbc->',
               g[v,v,o,o], g[o,v,o,o], g[o,o,v,o], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (pow(2,-1)*einsum('abij,ickl,jdac,klbd,ijab,jklacb,klbd->',
               g[v,v,o,o], g[o,v,o,o], g[o,v,v,v], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (pow(2,-1)*einsum('abij,ickl,kdab,jlcd,ijab,jklabc,jlcd->',
               g[v,v,o,o], g[o,v,o,o], g[o,v,v,v], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (pow(2,-1)*einsum('abij,cdak,ijcl,klbd,ijab,ijkcbd,klbd->',
               g[v,v,o,o], g[v,v,v,o], g[o,o,v,o], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (pow(2,-1)*einsum('abij,cdak,ikbl,jlcd,ijab,ijkbcd,jlcd->',
               g[v,v,o,o], g[v,v,v,o], g[o,o,v,o], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (einsum('abij,cdak,iebc,jkde,ijab,ijkbcd,jkde->',
               g[v,v,o,o], g[v,v,v,o], g[o,v,v,v], g[o,o,v,v], deltas['4'], deltas['6'], deltas['4'],
               optimize=True))
    mp += (pow(2,-1)*einsum('abij,ijkl,kcam,lmbc,ijab,klab,lmbc->',
               g[v,v,o,o], g[o,o,o,o], g[o,v,v,o], g[o,o,v,v], deltas['4'], deltas['4'], deltas['4'],
               optimize=True))
    mp += (pow(2,-1)*einsum('abij,icak,jklm,lmbc,ijab,jkbc,lmbc->',
               g[v,v,o,o], g[o,v,v,o], g[o,o,o,o], g[o,o,v,v], deltas['4'], deltas['4'], deltas['4'],
               optimize=True))
    mp += (-einsum('abij,icak,jdcl,klbd,ijab,jkcb,klbd->',
               g[v,v,o,o], g[o,v,v,o], g[o,v,v,o], g[o,o,v,v], deltas['4'], deltas['4'], deltas['4'],
               optimize=True))
    mp += (-einsum('abij,icak,kdbl,jlcd,ijab,jkbc,jlcd->',
               g[v,v,o,o], g[o,v,v,o], g[o,v,v,o], g[o,o,v,v], deltas['4'], deltas['4'], deltas['4'],
               optimize=True))
    mp += (pow(2,-1)*einsum('abij,icak,debc,jkde,ijab,jkbc,jkde->',
               g[v,v,o,o], g[o,v,v,o], g[v,v,v,v], g[o,o,v,v], deltas['4'], deltas['4'], deltas['4'],
               optimize=True))
    mp += (pow(2,-1)*einsum('abij,cdab,ieck,jkde,ijab,ijcd,jkde->',
               g[v,v,o,o], g[v,v,v,v], g[o,v,v,o], g[o,o,v,v], deltas['4'], deltas['4'], deltas['4'],
               optimize=True))
    return mp