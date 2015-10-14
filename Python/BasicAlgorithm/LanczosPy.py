from numpy import *
from scipy.sparse import csr_matrix
from scipy.linalg import eigh
from numpy.linalg import norm
from GlobalPy import RZERO,S_dtype
class Lanczos:
    '''
    The Lanczos class provides the following methods to deal with sparse Hermitian matrices:
    1) __init__: initialize an instance. Note that if the initial state to begin the iteration with is assigned, it must be normalized. A random/symmetric normalized state will be used as the initial state if it is not assigned.
    2) iter: the Lanczos iteration.
    3) tridiagnoal: returns the tridiagnoal matrix representation of the original sparse Hermitian matrix.
    4) eig: returns the ground state energy and optionally the ground state of the sparse Hermitian matrix.
    '''
    def __init__(self,matrix,vector=[],vtype='rd'):
        self.matrix=matrix
        if len(vector)==0:
            if vtype.lower()=='rd':
                self.new=zeros(matrix.shape[0],dtype=S_dtype)
                self.new[:]=random.rand(matrix.shape[0])
            else:
                self.new=ones(matrix.shape[0],dtype=S_dtype)
            self.new[:]=self.new[:]/norm(self.new)
        else:
            self.new=vector
        self.old=copy(self.new)
        self.cut=False
        self.a=[]
        self.b=[]

    def iter(self):
        count=len(self.a)
        buff=self.matrix.dot(self.new)
        self.a.append(vdot(self.new,buff))
        if count>0:
            buff[:]=buff[:]-self.a[count]*self.new-self.b[count-1]*self.old
        else:
            buff[:]=buff[:]-self.a[count]*self.new
        nbuff=norm(buff)
        if nbuff>RZERO:
            self.b.append(nbuff)
            self.old[:]=self.new[:]
            self.new[:]=buff[:]/nbuff
        else:
            self.cut=True
            self.b.append(0.0)
            self.old[:]=self.new[:]
            self.new[:]=0.0

    def tridiagnoal(self):
        nmatrix=len(self.a)
        result=zeros((nmatrix,nmatrix))
        for i,(a,b) in enumerate(zip(self.a,self.b)):
            result[i,i]=a.real
            if i<nmatrix-1: 
                result[i+1,i]=b
                result[i,i+1]=b
        return result

    def eig(self,job='n',precision=10**-10):
        if job in ('V','v'):gs=copy(self.new)
        delta=1.0;buff=inf
        while not self.cut and delta>precision:
            self.iter()
            if job in ('V','v'):
                w,vs=eigh(self.tridiagnoal())
                gse=w[0];v=vs[:,0]
            else:
                gse=eigh(self.tridiagnoal(),eigvals_only=True)[0]
            delta=abs(gse-buff)
            buff=gse
        if job in ('V','v'):
            self.a=[];self.b=[]
            for i in xrange(len(v)):
                if i==0:
                    self.new[:]=gs[:]
                    gs[:]=0.0
                gs[:]+=self.new*v[i]
                self.iter()
            return gse,gs
        else:
            return gse

# The following codes are used for tests only.
from scipy.sparse.linalg import eigsh
def test_lanczos():
    a=Lanczos(csr_matrix(array([[0.,-1.,-1.,0.],[-1.,0.,0.,-1.],[-1.,0.,0.,-1.],[0.,-1.,-1.,0.]])))
    print a.matrix.todense()
    print a.eig(job='v')
    print eigsh(a.matrix,which='SA',k=1)