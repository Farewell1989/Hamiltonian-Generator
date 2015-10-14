from Hamiltonian.TBAPy import *
from scipy.linalg import expm2,eig
class FLQT(TBA):
    '''
    Class FLQT deals with floquet problems.  
    '''
    def __init__(self,name=None,filling=0,mu=0,generator=None,**karg):
        super(FLQT,self).__init__(
            name=       name,
            filling=    filling,
            mu=         mu,
            generator=  generator
            )

    def evolution(self,t=[],**karg):
        nmatrix=len(self.generator.table)
        result=eye(nmatrix,dtype=GP.Q_dtype)
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
        savetxt(engine.dout+'/'+engine.name.full_name+'_EB.dat',result)
    if app.plot:
        plt.title(engine.name.full_name+'_EB')
        plt.plot(result[:,0],result[:,1:])
        if app.show:
            plt.show()
        else:
            plt.savefig(engine.dout+'/'+engine.name.full_name+'_EB.png')

# The following codes are used for tests only.
from BasicClass.BaseSpacePy import *
def test_flqt():
    N=50
    mu1=3.0
    mu2=0.0
    name='FLQT'
    p1=Point(site=0,rcoord=[0.0],icoord=[0.0],norbital=1,nspin=1,nnambu=2,scope=name)
    a1=array([1.0])
    a=FLQT(
        dout=       'Results',
        name=       name,
        generator=  Generator(
            lattice=    Lattice(name=name,points=[p1],translations=[(a1,N)]),
            #lattice=    Lattice(name=name,points=[p1],vectors=[a1]),
            parameters= {'mu1':mu1,'mu2':mu2},
            terms=[     Hopping('t1',-1.0),
                        Onsite('mu',0.0,modulate=lambda **karg: mu1 if karg['t']<1 else mu2),
                        Pairing('delta',0.5,neighbour=1,amplitude=lambda bond: 1 if bond.rcoord[0]>0 else -1)
                        ],
            nambu=True
            )
        )
    #a.addapps('EB',EB(path=BaseSpace('t',array([0,1])),save_data=False,run=TBAEB))
    a.addapps('EB',EB(ts=TSpace(array([0,1,2])),save_data=False,run=FLQTEB))
    a.runapps()
