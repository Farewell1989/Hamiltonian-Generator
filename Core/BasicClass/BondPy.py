'''
Bond.
'''
from BasicGeometryPy import *
class Bond:
    '''
    This class describes the bond in a lattice.
    Attributes:
        neighbour: integer
            The rank of the neighbour of the bond.
        spoint: Point
            the start point of the bond.
        epoint: Point
            the end point of the bond.
    '''
    
    def __init__(self,neighbour,spoint,epoint):
        '''
        Constructor.
        '''
        self.neighbour=neighbour
        self.spoint=spoint
        self.epoint=epoint

    def __str__(self):
        '''
        Convert an instance to string.
        '''
        return 'Neighbour: '+str(self.neighbour)+'\n'+'Start point:\n'+str(self.spoint)+'End point:\n'+str(self.epoint)
    
    @property
    def rcoord(self):
        '''
        The real coordinate of a bond.
        '''
        return self.epoint.rcoord-self.spoint.rcoord
    
    @property
    def icoord(self):
        '''
        The lattice coordinate of a bond.
        '''
        return self.epoint.icoord-self.spoint.icoord
    
    def is_intra_cell(self):
        '''
        Judge whether a bond is intra the unit cell or not. 
        '''
        if norm(self.icoord)< RZERO:
            return True
        else:
            return False

    @property
    def reversed(self):
        '''
        Return the reversed bond.
        '''
        return Bond(self.neighbour,self.epoint,self.spoint)
