from IndexPy import *
from OperatorPy import *
from TermPy import *
from BondPy import *
from TablePy import *
import GlobalPy as GP
class Hubbard(Term):
    '''
    Class Hubbard provides a complete description of single orbital and multi orbital Hubbard interactions.
    1) value: an array with len=1 for single orbital systems and len=3 for multi orbital systems. For multi orbital systems, the interactions include intra-orbital, inter-orbital, spin flip and pair hopping terms.
    2) atom: the atom index of the point where the Hubbard interactions are defined.
    '''
    def __init__(self,tag,value,atom=None,modulate=None):
        super(Hubbard,self).__init__('hb',tag,value,modulate)
        if not atom is None: self.atom=atom

    def __str__(self):
        '''
        Convert an instance to string.
        '''
        result=''
        if hasattr(self,'atom'): result+='Atom: '+str(self.atom)+'\n'
        result+='Tag,value: '+self.tag+','+str(self.value)+'\n'
        return result

    def __add__(self,other):
        '''
        Overloaded operator(+), which supports the addition of a Hubbard instance with a Hubbard/HubbardList instance.
        '''
        result=HubbardList(deepcopy(self))
        if isinstance(other,Hubbard):
            result.append(deepcopy(other))
        elif isinstance(other,HubbardList):
            result.extend(deepcopy(other))
        else:
            raise ValueError('Hubbard "+" error: the other parameter must be an instance of Hubbard or HubbardList')
        return result

    def __pos__(self):
        '''
        Overloaded operator(+), i.e. +self.
        '''
        result=HubbardList(deepcopy(self))
        return result

    def mesh(self,bond):
        '''
        Generate the mesh of Hubbard terms.
        '''
        ndim=bond.epoint.norbital*bond.epoint.nspin
        result=zeros((ndim,ndim,ndim,ndim),dtype=GP.H_dtype)
        if hasattr(self,'atom'):
            atom=self.atom
        else:
            atom=bond.epoint.atom
        if atom==bond.epoint.atom:
            nv=len(self.value)
            if nv>=1:
                for h in xrange(bond.epoint.norbital):
                    i=bond.epoint.seq_state(h,1,ANNIHILATION)
                    j=bond.epoint.seq_state(h,0,ANNIHILATION)
                    k=j
                    l=i
                    result[i,j,k,l]=self.value[0]/2
            if nv==3:
                for h in xrange(bond.epoint.norbital):
                    for g in xrange(bond.epoint.norbital):
                      if g!=h:
                        i=bond.epoint.seq_state(g,1,ANNIHILATION)
                        j=bond.epoint.seq_state(h,0,ANNIHILATION)
                        k=j
                        l=i
                        result[i,j,k,l]=self.value[1]/2
                for h in xrange(bond.epoint.norbital):
                    for g in xrange(h):
                        for f in xrange(2):
                            i=bond.epoint.seq_state(g,f,ANNIHILATION)
                            j=bond.epoint.seq_state(h,f,ANNIHILATION)
                            k=j
                            l=i
                            result[i,j,k,l]=(self.value[1]-self.value[2])/2
                for h in xrange(bond.epoint.norbital):
                    for g in xrange(h):
                        i=bond.epoint.seq_state(g,1,ANNIHILATION)
                        j=bond.epoint.seq_state(h,0,ANNIHILATION)
                        k=bond.epoint.seq_state(g,0,ANNIHILATION)
                        l=bond.epoint.seq_state(h,1,ANNIHILATION)
                        result[i,j,k,l]=self.value[2]
                for h in xrange(bond.epoint.norbital):
                    for g in xrange(h):
                        i=bond.epoint.seq_state(g,1,ANNIHILATION)
                        j=bond.epoint.seq_state(g,0,ANNIHILATION)
                        k=bond.epoint.seq_state(h,0,ANNIHILATION)
                        l=bond.epoint.seq_state(h,1,ANNIHILATION)
                        result[i,j,k,l]=self.value[2]
        return result

class HubbardList(list):
    '''
    The HubbardList class pack several Hubbard instances as a whole for convenience.
    '''
    
    def __init__(self,*arg):
        self.mode='hb'
        for obj in arg:
            if isinstance(obj,Hubbard):
                self.append(obj)
            else:
                raise ValueError("HubbardList init error: the input argument should be Hubbard instances.")

    def __str__(self):
        '''
        Convert an instance to string.
        '''
        result='Hubbard terms:\n'
        for i,v in enumerate(self):
            result+='Node '+str(i)+':\n'+str(v)

    def __add__(self,other):
        '''
        Overloaded operator(+), which supports the addition of a HubbardList instance with a Hubbard/HubbardList instance.
        '''
        result=HubbardList(*deepcopy(self))
        if isinstance(other,Hubbard):
            result.append(deepcopy(other))
        elif isinstance(other,HubbardList):
            result.extend(deepcopy(other))
        else:
            raise ValueError('HubbardList "+" error: the other parameter must be an instance of Hubbard or HubbardList')
        return result

    def mesh(self,bond):
        '''
        Generate the mesh of all Hubbard terms defined on a bond.
        '''
        if bond.epoint.nspin==2 and bond.spoint.nspin==2:
            ndim=bond.epoint.norbital*bond.epoint.nspin
            result=zeros((ndim,ndim,ndim,ndim),dtype=GP.H_dtype)
            if bond.neighbour==0:
                for obj in self:
                    result+=obj.mesh(bond)
            return result
        else:
            raise ValueError('HubbardList mesh error: the input bond must be onsite ones nspin=2.')

    def operators(self,bond,table,half=True):
        '''
        Generate the set of non-zero operators defined on the input bond.
        1) The index sequences are determined by the index-sequence table. Since Hubbard terms are quartic and cannot be represented in lattice representations, the sequences will never be in the Nambu space.
        2) Because of the hermiticity of the Hamiltonian, when the parameter 'half' is set to be true, only one half of the set is returned. Note that the coefficient of the self-hermitian operators is also divided by a factor 2 so that the whole set exactly equals the returned set plus its Hermitian conjugate. The half=False case is not supported yet.
        '''
        result=OperatorList()
        buff=self.mesh(bond)
        indices=argwhere(abs(buff)>RZERO)
        for (i,j,k,l) in indices:
            index1=Index(scope=bond.epoint.scope,site=bond.epoint.site,**bond.epoint.state_index(i))
            index2=Index(scope=bond.epoint.scope,site=bond.epoint.site,**bond.epoint.state_index(j))
            index3=Index(scope=bond.epoint.scope,site=bond.epoint.site,**bond.epoint.state_index(k))
            index4=Index(scope=bond.epoint.scope,site=bond.epoint.site,**bond.epoint.state_index(l))
            result.append(E_Hubbard(buff[i,j,k,l],indices=deepcopy([index1.dagger,index2.dagger,index3,index4]),rcoords=[bond.epoint.rcoord],icoords=bond.epoint.icoord,seqs=[table[index1],table[index2],table[index3],table[index4]]))
        return result

# The following codes are used for tests only.
from LatticePy import *
def test_hubbard():
    p1=Point(site=0,rcoord=[0.0,0.0],icoord=[0,0],norbital=2,nspin=2,nnambu=1,scope="WG")
    l=Lattice(name="WG",points=[p1])
    l.plot(show='y')
    table=Table(l.indices(nambu=True))
    a=HubbardList(Hubbard('U,U,J',[20.0,12.0,5.0]))
    opts=OperatorList()
    for bonds in l.bonds:
        for bond in bonds:
            opts.extend(a.operators(bond,table))
    print opts