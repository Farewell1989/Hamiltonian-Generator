from Hamiltonian.Core.BasicClass import *
class TriangleDataBase:
    def __init__(self,scope,norbital=1,nspin=2,nnambu=1):
        self.points=[]
        self.vectors=[]
        if scope=='T1':
            self.points.append(Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],norbital=norbital,nspin=nspin,nnambu=nnambu,scope=scope))
            self.vectors.append(array([1.0,0.0]))
            self.vectors.append(array([0.5,sqrt(3)/2]))
        elif scope=='T12':
            self.points.append(Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],norbital=norbital,nspin=nspin,nnambu=nnambu,scope=scope))
            self.points.append(Point(site=1,rcoord=[1.0,0.0],icoord=[0.0,0.0],norbital=norbital,nspin=nspin,nnambu=nnambu,scope=scope))
            self.points.append(Point(site=2,rcoord=[2.0,0.0],icoord=[0.0,0.0],norbital=norbital,nspin=nspin,nnambu=nnambu,scope=scope))
            self.points.append(Point(site=3,rcoord=[3.0,0.0],icoord=[0.0,0.0],norbital=norbital,nspin=nspin,nnambu=nnambu,scope=scope))
            self.points.append(Point(site=4,rcoord=[0.5,-sqrt(3.0)/2],icoord=[0.0,0.0],norbital=norbital,nspin=nspin,nnambu=nnambu,scope=scope))
            self.points.append(Point(site=5,rcoord=[1.5,-sqrt(3.0)/2],icoord=[0.0,0.0],norbital=norbital,nspin=nspin,nnambu=nnambu,scope=scope))
            self.points.append(Point(site=6,rcoord=[2.5,-sqrt(3.0)/2],icoord=[0.0,0.0],norbital=norbital,nspin=nspin,nnambu=nnambu,scope=scope))
            self.points.append(Point(site=7,rcoord=[1.0,-sqrt(3.0)],icoord=[0.0,0.0],norbital=norbital,nspin=nspin,nnambu=nnambu,scope=scope))
            self.points.append(Point(site=8,rcoord=[2.0,-sqrt(3.0)],icoord=[0.0,0.0],norbital=norbital,nspin=nspin,nnambu=nnambu,scope=scope))
            self.points.append(Point(site=9,rcoord=[0.5,sqrt(3.0)/2],icoord=[0.0,0.0],norbital=norbital,nspin=nspin,nnambu=nnambu,scope=scope))
            self.points.append(Point(site=10,rcoord=[1.5,sqrt(3.0)/2],icoord=[0.0,0.0],norbital=norbital,nspin=nspin,nnambu=nnambu,scope=scope))
            self.points.append(Point(site=11,rcoord=[2.5,sqrt(3.0)/2],icoord=[0.0,0.0],norbital=norbital,nspin=nspin,nnambu=nnambu,scope=scope))
            self.vectors.append(array([0.0,2*sqrt(3.0)]))
            self.vectors.append(array([3.0,sqrt(3.0)]))
        else:
            raise ValueError('TriangleDataBase error: the scope parameter not supported.')

if __name__=='__main__':
      for scope in ['T1','T12']:
        buff=TriangleDataBase(scope)
        l=Lattice(name=scope,points=buff.points,vectors=buff.vectors)
        l.plot()
