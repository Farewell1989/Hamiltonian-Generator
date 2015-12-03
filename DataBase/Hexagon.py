from Hamiltonian.Core.BasicClass import *
class HexagonDataBase:
    def __init__(self,scope,norbital=1,nspin=2,nnambu=1):
        A=Fermi(atom=0,norbital=norbital,nspin=nspin,nnambu=nnambu)
        B=Fermi(atom=1,norbital=norbital,nspin=nspin,nnambu=nnambu)
        self.points=[]
        self.vectors=[]
        if scope=='H2':
            self.points.append(Point(scope=scope,site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=1,rcoord=[0.0,sqrt(3)/3],icoord=[0.0,0.0],struct=B))
            self.vectors.append(array([1.0,0.0]))
            self.vectors.append(array([0.5,sqrt(3)/2]))
        elif scope=='H4':
            self.points.append(Point(scope=scope,site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=1,rcoord=[0.0,sqrt(3)/3],icoord=[0.0,0.0],struct=B))
            self.points.append(Point(scope=scope,site=2,rcoord=[0.5,sqrt(3)/2],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=3,rcoord=[0.5,-sqrt(3)/6],icoord=[0.0,0.0],struct=B))
            self.vectors.append(array([1.0,0.0]))
            self.vectors.append(array([0.0,sqrt(3.0)]))        
        elif scope=='H6':
            self.points.append(Point(scope=scope,site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=1,rcoord=[0.0,sqrt(3)/3],icoord=[0.0,0.0],struct=B))
            self.points.append(Point(scope=scope,site=2,rcoord=[0.5,sqrt(3)/2],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=3,rcoord=[0.5,-sqrt(3)/6],icoord=[0.0,0.0],struct=B))
            self.points.append(Point(scope=scope,site=4,rcoord=[1.0,0.0],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=5,rcoord=[1.0,sqrt(3)/3],icoord=[0.0,0.0],struct=B))
            self.vectors.append(array([1.5,sqrt(3)/2]))
            self.vectors.append(array([1.5,-sqrt(3)/2]))
        elif scope=='H8P':
            self.points.append(Point(scope=scope,site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=1,rcoord=[0.0,sqrt(3)/3],icoord=[0.0,0.0],struct=B))
            self.points.append(Point(scope=scope,site=2,rcoord=[0.5,sqrt(3)/2],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=3,rcoord=[0.5,-sqrt(3)/6],icoord=[0.0,0.0],struct=B))
            self.points.append(Point(scope=scope,site=4,rcoord=[1.0,0.0],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=5,rcoord=[1.0,sqrt(3)/3],icoord=[0.0,0.0],struct=B))
            self.points.append(Point(scope=scope,site=6,rcoord=[0.5,-sqrt(3)/2],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=7,rcoord=[0.5,sqrt(3)*5/6],icoord=[0.0,0.0],struct=B))
            self.vectors.append(array([1.0,sqrt(3)]))
            self.vectors.append(array([1.5,-sqrt(3)/2]))
        elif scope=='H8O':
            self.points.append(Point(scope=scope,site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=1,rcoord=[0.0,sqrt(3)/3],icoord=[0.0,0.0],struct=B))
            self.points.append(Point(scope=scope,site=2,rcoord=[0.5,sqrt(3)/2],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=3,rcoord=[0.5,-sqrt(3)/6],icoord=[0.0,0.0],struct=B))
            self.points.append(Point(scope=scope,site=4,rcoord=[1.0,0.0],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=5,rcoord=[1.0,sqrt(3)/3],icoord=[0.0,0.0],struct=B))
            self.points.append(Point(scope=scope,site=6,rcoord=[1.5,sqrt(3)/2],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=7,rcoord=[1.5,-sqrt(3)/6],icoord=[0.0,0.0],struct=B))
            self.vectors.append(array([2.0,0.0]))
            self.vectors.append(array([0.0,sqrt(3)]))
        elif scope=='H10':
            self.points.append(Point(scope=scope,site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=1,rcoord=[0.0,sqrt(3)/3],icoord=[0.0,0.0],struct=B))
            self.points.append(Point(scope=scope,site=2,rcoord=[0.5,sqrt(3)/2],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=3,rcoord=[0.5,-sqrt(3)/6],icoord=[0.0,0.0],struct=B))
            self.points.append(Point(scope=scope,site=4,rcoord=[1.0,0.0],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=5,rcoord=[1.0,sqrt(3)/3],icoord=[0.0,0.0],struct=B))
            self.points.append(Point(scope=scope,site=6,rcoord=[1.5,sqrt(3)/2],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=7,rcoord=[1.5,-sqrt(3)/6],icoord=[0.0,0.0],struct=B))
            self.points.append(Point(scope=scope,site=8,rcoord=[2.0,0.0],icoord=[0.0,0.0],struct=A))
            self.points.append(Point(scope=scope,site=9,rcoord=[2.0,sqrt(3)/3],icoord=[0.0,0.0],struct=B))
            self.vectors.append(array([2.5,sqrt(3)/2]))
            self.vectors.append(array([0.0,sqrt(3)]))
        else:
            raise ValueError('HexagonDataBase error: the scope parameter not supported.')

if __name__=='__main__':
    for scope in ['H2','H4','H6','H8P','H8O','H10']:
        buff=HexagonDataBase(scope)
        l=Lattice(name=scope,points=buff.points,vectors=buff.vectors)
        l.plot()
