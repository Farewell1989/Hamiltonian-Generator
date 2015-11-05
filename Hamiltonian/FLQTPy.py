from Hamiltonian.TBAPy import *
from scipy.linalg import expm2,eig
class FLQT(TBA):
    '''
    Class FLQT deals with floquet problems.  
    '''
    def __init__(self,filling=0,mu=0,lattice=None,terms=None,nambu=False,**karg):
        super(FLQT,self).__init__(
            filling=    filling,
            mu=         mu,
            lattice=    lattice,
            terms=      terms,
            nambu=      nambu
            )

    def evolution(self,t=[],**karg):
        nmatrix=len(self.generator.table)
        result=eye(nmatrix,dtype=complex128)
        nt=len(t)
        for i,time in enumerate(t):
            if i<nt-1:
                result=dot(expm2(-1j*self.matrix(t=time,**karg)*(t[i+1]-time)),result)
        return result

def FLQTEB(engine,app):
    nmatrix=len(engine.generator.table)
    if app.path!=None:
        result=zeros((app.path.rank,nmatrix+1))
        if len(app.path.mesh.shape)==1:
            result[:,0]=app.path.mesh
        else:
            result[:,0]=array(xrange(app.path.rank))
        for i,parameter in enumerate(list(app.path.mesh)):
            result[i,1:]=phase(eig(engine.evolution(t=app.ts.mesh,**{app.path.mode:parameter}))[0])/app.ts.volume
    else:
        result=zeros((2,nmatrix+1))
        result[:,0]=array(xrange(2))
        result[0,1:]=angle(eig(engine.evolution(t=app.ts.mesh))[0])/app.ts.volume
        result[1,1:]=result[0,1:]
    if app.save_data:
        savetxt(engine.dout+'/'+engine.name.full+'_EB.dat',result)
    if app.plot:
        plt.title(engine.name.full+'_EB')
        plt.plot(result[:,0],result[:,1:])
        if app.show:
            plt.show()
        else:
            plt.savefig(engine.dout+'/'+engine.name.full+'_EB.png')
