from __future__ import division
import numpy as np

if __name__ == '__main__':

    '''
    import the HUGENHOLTZ class from source.hugenholtz and run for the required
    order as
    '''
    from source.hugenholtz import HUGENHOLTZ

    order = 4
    h = HUGENHOLTZ(order)

    '''
    the number of node pairs, number of lines and total number of diagrams can
    be obtained by
    '''
    print('Order {:2} has {:3} node pairs, {:3} lines between the nodes, and a total of {:3} diagrams'.
          format(h.order, h.count_nodal_pairs(), h.count_nodal_lines(), h.diagrams))

    '''
    a list of the nodal pairs can be obtained by - self.nodal_pairs
    '''

    '''
    the list of particle (up arrow) connections and hole (down arrow) connections can be obtained
    by... each connection is itself a list each element being the number of up/down connections
    between the corresponding node pair
    '''
    holes = h.holes()
    print('\n diagram    up connections           down connections\n-------------------------------------------------------')
    for diagram in range(h.diagrams):
        print('   ', str(diagram+1).ljust(3,' '), '  ', h.particles[diagram][1],'     ', holes[diagram][1])

    '''
    The rules for interpreting the Hugenholtz diagrams are
    1. Each node contributes an antisymmetrized integral of the form < || > to the numerator. The elements of the integral
       are < in_1 in_2 || out_1 out_2 >. (order is arbitary and accounted for by the phase calculation)

    2. Each adjacent pair of nodes contributes to the denominator a factor (sum of hole lines) - (sum of particle lines). The
       sums run over all lines cutting an imaginary horizontal line cutting the vertical between the two adjacent nodes.

    3. If h is the number of down lines and l is the number of closed loops then the overall sign is (-1)**(h+l). Closed loops
       are determined by writing out the indices in order from 1. eg <ab||rs> <rs||ab>, then start at first index 'a' and skip
       along the list until you get back to an 'a' (the index list wraps): sp 'a'->'r'->'r'->'a'. Now go to the next label we
       haven't yet visited 'b', and repeat the skipping operation 'b'->'s'->'s'->'b'. We have now visited every index so the
       closed loops are 2 - the number of skipping operations completed.

    4. There is a weight factor which the expression is to be multiplied by which is (2)**(-k), where k are the number of
       pairs of equivalent lines. A pair of lines are equivalent if they begin and end at the same apir of nodes in the same
       direction.
       
    5. Sum the diagram over each index

    How these rules are computed can be seen for a diagram by - to find the number of a diagram either look at the .tex file
    or a plot on which the numbers are given. Alternatively if the up and down connections are known the number can be found
    from h.ordinal(([0,1,1,1,1,0], [0,0,2,2,0,0])) [h,ordinal(up, down)).
    '''
    diagram = 1
    h.debug(1)

    '''
    we see from this output that how the numerator and denominator are calculated together with the factors
    '''

    '''
    to create a .tex file of the algrebraic expressions use
    '''
    h.generate_file('latex', file='my_algebra_mp4.tex')

    '''
    to create a .py file which can be used to evaluate the MPn energy correction. The subroutine produced will require
    the orbital energies and 2-electron repulsion integrals in a molecular spin-basis, and the number of occupied spin
    orbitals
    '''
    h.generate_file('python', file='my_python.py')

    '''
    it is possible to get a plot of the diagrams from the DISPLAY class. An individual diagram can be plotted by
    '''
    from source.display import DISPLAY
    DISPLAY(h, diagram=1, ident=True)

    '''
    if ident is True (default False) then the op and down connection will be plotted on the diagram. The diagram number is always 
    plotted. Ident can only be enabled for individual diagrams. If the diagram is set to None (which is default) all diagrams for 
    the specified order will be plotted. Group plotting is only available for orders up to and including 4 as beyond that the 
    number becomes unmanagable. Below is an example of plotting all diagrams for specified order and save the resulting plot.
    '''
    DISPLAY(h, diagram=None, file='my_plot_mp4.png')
    
