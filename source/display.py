from __future__ import division
import matplotlib.pyplot as mp
from matplotlib.patches import Arc
import numpy as np

class DISPLAY(object):
    '''
    class to draw a HUGENHOLTZ object as a directed graph
    '''

    def __init__(self, HUGENHOLTZ, diagram=None, ident=False, file=None):

        self.h = HUGENHOLTZ
        self.diagram, self.ident = diagram, ident

        #special case order 2
        if self.h.order == 2: self.diagram = 1
        
        size = 3
        if type(diagram) == tuple:
            self.diagram = self.h.ordinal(self.diagram)
        if not self.diagram is None:
            fig, ax = mp.subplots(figsize=(size, size))
            fig.suptitle('Connected Hugenholtz diagrams for order ' + str(self.h.order), fontsize='small')
            self.draw_diagram(ax) 
        else:
            #for group plot reduce size and switch off idents
            size, self.ident = 2, False
            if self.h.order == 2:
                ncol, nrow, height, width = 1, 1, 3, 3
            if self.h.order == 3:
                ncol, nrow, height, width = 3, 1, 3, 5
            if self.h.order == 4:
                ncol, nrow, height, width = 6, 7, 12, 10
            if self.h.order > 4:
                print('individual diagrams only for order > 4')
                return
            
            fig, axs = mp.subplots(nrow, ncol, figsize=(width, height))
            fig.suptitle('Connected Hugenholtz diagrams for order ' + str(self.h.order))
            for i, ax in enumerate(axs.ravel()):
                self.diagram = i + 1
                if i >= self.h.diagrams:
                    ax.axis('off')
                    ax.plot()
                else:
                    self.draw_diagram(ax)

        #savefile if requested
        if not file is None:
            mp.savefig(file)
                    
        mp.show()

        return
                              
    def draw_diagram(self, ax):
        '''
        draw directed graph of 'diagram'
        '''

        #get diagram statistics
        up, down = self.h.particles[self.diagram - 1][1], self.h.holes()[self.diagram - 1][1]
                              
        #frame size and internodal distance with top and bottom margins 0.5 internode distance
        frame_size = 10
        vertical_separation = frame_size/self.h.order
        vertical_margin = vertical_separation * 0.5

        #node coordinates
        node_x = [0] * self.h.order
        node_y = [i*vertical_separation+vertical_margin for i in range(self.h.order)]

        #initialise plot - center line of diagram at x=0
        ax.set_xlim([-frame_size//2, frame_size//2])
        ax.set_ylim([0, frame_size])

        if self.ident:
            ax.text(-5, frame_size - vertical_margin, '[' + ','.join([str(x) for x in up])   + '] \u2191', fontsize='small')
            ax.text(-5,              vertical_margin, '[' + ','.join([str(x) for x in down]) + '] \u2193', fontsize='small', va='top')
        ax.text(frame_size - 5.0, frame_size, str(self.diagram), fontsize='xx-small', style='italic',
                                              bbox={'facecolor': 'orange', 'alpha': 0.5, 'boxstyle':'Circle,pad=0.2'})
        
        #plot node points
        ax.plot(node_x, node_y, '.k')
        ax.axis('off')

        #marker size
        marker_size = 3.5 if self.h.order >= 4 else 5
        
        def ellipsoid_arc(center, radius, side_left=True, arrow_up=True, e=1):
            '''
            plot a elliptical arc center at (x, y), width and height (w, h) and angle subtended (theta)
            '''
            x, y = center
    
            #size of ellispe
            major_axis = (frame_size/self.h.order) * radius
            minor_axis = major_axis * 0.75 * e

            #subtended angle
            theta = [90, 270] if side_left else [-90, 90]

            #arrow position
            arrow_y = y
            arrow_x = x - minor_axis * 0.5 if side_left else x + minor_axis * 0.5
            pointer = ['^', 'r'] if arrow_up else ['v', 'b']
            ax.plot(arrow_x, arrow_y, marker=pointer[0], color=pointer[1], markersize=marker_size)

            h = Arc((x, y), minor_axis, major_axis, theta1= theta[0], theta2 = theta[1], linewidth=1)

            return h

        #pairs - tuple of node end-points for each pair
        pairs = [(i,j) for i in range(self.h.order-1) for j in range(i+1, self.h.order)]

        #pair properties - _flow (up, down) line count between each pair, _order is number of lines between each node pair
        #_state - the initial state (up lines, down lines, draw left(True) draw right(False)
        pairs_flow = [(a, b) for a, b in zip(up, down)]             
        pairs_order = [sum(i) for i in pairs_flow]
        pairs_state = [[j, j, True] for i in range(self.h.order-1, -1, -1) for j in range(i) ]

        #node number at which new  node sequence starts
        base_nodes = [0] 
        x = 0
        for i in range(self.h.order-1, 1, -1):
            x += i
            base_nodes.append(x)

        #current node
        node = 0                                                                     
        asymmetric_direction = True

        while True:

            #values for current node
            f, t = pairs[node][0], pairs[node][1]
            y = [vertical_margin + f * vertical_separation, vertical_margin + t * vertical_separation]
            center = (y[0] + y[1]) * 0.5
    
            if pairs_order[node] == 0:
                pass
    
            if pairs_order[node] == 1:
                #if base pair no obstruction draw straight line between nodes
                if node in base_nodes:
                    ax.plot([0.0, 0.0], y, linewidth=2, color='k')
                    #up arrow if flow (*, ) is 1
                    pointer = ['^', 'r'] if pairs_flow[node][0] == 1 else ['v', 'b']
                    ax.plot(0.0, center, marker=pointer[0], color=pointer[1], markersize=marker_size)
            
                else:
                    #there is another node between this node pair - draw arc around it
                    if pairs_state[node][0] != pairs_state[node][1]:
                        left = True if pair_state[node][0] <= pairs_state[node][1] else False
                    else:
                        left = True if asymmetric_direction else False
                
                    up = (pairs_flow[node][0] != 0)
                
                    #draw arc round intervening nodes and update state information for this pair
                    ax.add_patch(ellipsoid_arc((0, center), pairs_state[node][0] + 1, side_left=left, arrow_up=up))                       

                    pairs_state[node][i] += 1
                    pairs_state[node][2] = not pairs_state[node][2]
                    asymmetric_diection = pairs_state[node][2]
                
            if pairs_order[node] == 2:
                #draw pair of arcs at appropriate radii - one up and one down on opposite sides
                up = (pairs_flow[node][0] != 0)
                ax.add_patch(ellipsoid_arc((0, center), pairs_state[node][0] + 1, side_left=True,  arrow_up=up))
                up = (pairs_flow[node][1] == 0)
                ax.add_patch(ellipsoid_arc((0, center), pairs_state[node][1] + 1, side_left=False, arrow_up=up))                       
       
                pairs_state[node][0] += 1
                pairs_state[node][1] += 1
  
            if pairs_order[node] == 3:
                
                #if first pair in sequence can drwa straight line between nodes and pair of arcs
                if node in base_nodes:

                    ax.plot([0, 0], y, linewidth=2, color='k')
                    pointer = ['^', 'r'] if pairs_flow[node][0] != 1 else ['v', 'b']
                    ax.plot(0, center, marker=pointer[0], color=pointer[1], markersize=marker_size)
                    
                    ax.add_patch(ellipsoid_arc((0, center), pairs_state[node][0] + 1, side_left=True,  arrow_up=True))
                    ax.add_patch(ellipsoid_arc((0, center), pairs_state[node][1] + 1, side_left=False, arrow_up=False))                       

                else:
                    #intervening node(s) draw pair and extra arc around intervening node(s)
                    ax.add_patch(ellipsoid_arc((0, center), pairs_state[node][0] + 1, side_left=True,  arrow_up=True))
                    ax.add_patch(ellipsoid_arc((0, center), pairs_state[node][1] + 1, side_left=False, arrow_up=False))                       

                    if pairs_state[node][0] != pairs_state[node][1]:
                        left = True if pair_state[node][0] <= pair_state[node][1] else False
                    else:
                        left = True if asymmetric_direction else False
                        
                    up = (pairs_flow[node][0] > pairs_flow[node][1])
                    ax.add_patch(ellipsoid_arc((0, center), pairs_state[node][0] + 1, side_left=left, arrow_up=up, e=0.5))

                    pairs_state[node][i] += 1
                    pairs_state[node][2] = not pairs_state[node][2]
                    asymmetric_direction = pairs_state[node][2]

            #special case for order 2
            if pairs_order[node] == 4:
                ax.add_patch(ellipsoid_arc((0, center), pairs_state[node][0] + 1, side_left=True,  arrow_up=True))
                ax.add_patch(ellipsoid_arc((0, center), pairs_state[node][1] + 1, side_left=False, arrow_up=False))                       
                ax.add_patch(ellipsoid_arc((0, center), pairs_state[node][0] + 1, side_left=True,  arrow_up=True, e=0.5))
                ax.add_patch(ellipsoid_arc((0, center), pairs_state[node][1] + 1, side_left=False, arrow_up=False, e=0.5))                       

            #go to nect node pairunless finished  
            node += 1
            asymmetric_direction = not asymmetric_direction

            if node == len(pairs_order): break

            
