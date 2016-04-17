# -*- coding: utf-8 -*-
"""
Created on Fri Jan 08 16:09:37 2016

@author: DZY
"""

'''
Linear Spin wave theory
'''
from Hamiltonian.Core.BasicClass.AppPackPy import *
from Hamiltonian.Core.BasicClass.NamePy import *
from Hamiltonian.Core.BasicAlgorithm.BerryCurvaturePy import *
from Hamiltonian.Core.BasicClass.DZYOrderPy import *
from Hamiltonian.Core.BasicClass.LatticePy import *
from Hamiltonian.Core.BasicClass.TablePy import *
from Hamiltonian.Core.BasicClass.DZYSuperexchangePy import *
from Hamiltonian.Core.BasicClass.DZYSpinWGeneratorPy import *
from Hamiltonian.Core.BasicClass.DZYSpinWmatrix import *
from Hamiltonian.Core.BasicClass.DZYSuperbondPy import *
from Hamiltonian.Core.BasicClass.DZYSindexPackPy import *
from scipy.linalg import eig
import matplotlib.pyplot as plt
from scipy.io import savemat  

class SpinW(Engine):
    
    
    def __init__(self,lattice=None,terms=None,order=None,**karg):
        '''
        Constructor.
        '''
        self.lattice=lattice
        self.terms=terms
        self.order=order
        self.generators={}
        self.generators['h']=SpinWGenerator(self.lattice.dim,self.order,self.terms.bondclassify(self.lattice.superbonds),terms=self.terms)
        self.name.update(const=self.generators['h'].parameters['const'])

    def matrix(self,k=[],**karg):
        '''
        This method returns the matrix representation of the Hamiltonian.
        Parameters:
            k: 1D array-like, optional
                The coords of a point in K-space.
            karg: dict, optional
                Other parameters.
        Returns:
            result: 2D ndarray
                The matrix representation of the Hamiltonian.
        '''
        self.generators['h'].update(**karg)
        m=len(self.generators['h'].table)
        result=zeros((2*m,2*m),dtype=complex128)
        for opt in self.generators['h'].operators:
            phase=1 if len(k)==0 else exp(-1j*inner(k,opt.rcoords[0]))
            result[opt.seqs]+=opt.value*phase
            i,j=opt.seqs
            result[(mod(j+m,2*m),mod(i+m,2*m))]+=opt.value*conjugate(phase)
        result+=conjugate(result.T)
        return kron(array([[1.0,0.0],[0.0,-1.0]]),eye(m)).dot(result)

    def matrices(self,basespace=None,mode='*'):
        '''
        This method returns a generator which iterates over all the Hamiltonians living on the input basespace.
        Parameters:
            basespace: BaseSpace,optional
                The base space on which the Hamiltonians lives.
            mode: string,optional
                The mode which the generators takes to iterate over the base space.
        Returns:
            yield a 2D ndarray.
        '''
        if basespace is None:
            yield self.matrix()
        else:
            for paras in basespace(mode):
                yield self.matrix(**paras)

    def eigvals(self,basespace=None):
        '''
        This method returns all the eigenvalues of the Hamiltonian.
        Parameters:
            basespace: BaseSpace, optional
                The base space on which the Hamiltonian is defined.
        Returns:
            result: 1D ndarray
                All the eigenvalues.
        '''
        nmatrix=2*len(self.generators['h'].table)
        result=zeros(nmatrix*(1 if basespace==None else product(basespace.rank.values())))
        if basespace is None:
            result[...]=eigh(self.matrix(),eigvals_only=True)
        else:
            for i,paras in enumerate(basespace('*')):
                result[i*nmatrix:(i+1)*nmatrix]=eigh(self.matrix(**paras),eigvals_only=True)
        return result


def TBAEB(engine,app):
    nmatrix=2*len(engine.generators['h'].table)
    if app.path!=None:
        key=app.path.mesh.keys()[0]
        result=zeros((app.path.rank[key],nmatrix+1),dtype=complex128)
        if len(app.path.mesh[key].shape)==1:
            result[:,0]=app.path.mesh[key]
        else:
            result[:,0]=array(xrange(app.path.rank[key]))
        for i,parameter in enumerate(list(app.path.mesh[key])):
            result[i,1:]=eig(engine.matrix(**{key:parameter}),right=False)
    else:
        result=zeros((2,nmatrix+1),dtype=complex128)
        result[:,0]=array(xrange(2))
        result[0,1:]=eig(engine.matrix(),right=False)
        result[1,1:]=result[0,1:]
    if app.save_data:
        savemat(engine.dout+'/'+'_EB.dat',{'out':result})
    if app.plot:
        plt.figure(1)
        plt.title(engine.name.full+'_EB')
        plt.plot(result[:,0],sort(real(result[:,1:])))
        plt.figure(2)
        plt.title(engine.name.full+'_IMAG')
        plt.plot(result[:,0],sort(imag(result[:,1:])))
        if app.show:
            plt.show()
        else:
            plt.savefig(engine.dout+'/'+engine.name.full+'_EB.png')
        plt.close()

#def TBADOS(engine,app):
#    result=zeros((app.ne,2))
#    eigvals=engine.eigvals(app.BZ)
#    for i,v in enumerate(linspace(eigvals.min(),eigvals.max(),num=app.ne)):
#       result[i,0]=v
#       result[i,1]=sum(app.eta/((v-eigvals)**2+app.eta**2))
#    if app.save_data:
#        savetxt(engine.dout+'/'+engine.name.full+'_DOS.dat',result)
#    if app.plot:
#        plt.title(engine.name.full+'_DOS')
#        plt.plot(result[:,0],result[:,1])
#        if app.show:
#            plt.show()
#        else:
#            plt.savefig(engine.dout+'/'+engine.name.full+'_DOS.png')
#        plt.close()

#def TBACP(engine,app):
#    nelectron=int(round(engine.filling*app.BZ.rank['k']*len(engine.generators['h'].table)))
#    eigvals=sort((engine.eigvals(app.BZ)))
#    app.mu=(eigvals[nelectron]+eigvals[nelectron-2])/2
#    engine.mu=app.mu
#    print 'mu:',app.mu
#
#def TBACN(engine,app):
#    H=lambda kx,ky: engine.matrix(k=[kx,ky])
#    app.bc=zeros(app.BZ.rank['k'])
#    for i,paras in enumerate(app.BZ()):
#        app.bc[i]=berry_curvature(H,paras['k'][0],paras['k'][1],engine.mu,d=app.d)
#    print 'Chern number(mu):',app.cn,'(',engine.mu,')'
#    if app.save_data or app.plot:
#        buff=zeros((app.BZ.rank['k'],3))
#        buff[:,0:2]=app.BZ.mesh['k']
#        buff[:,2]=app.bc
#    if app.save_data:
#        savetxt(engine.dout+'/'+engine.name.full+'_BC.dat',buff)
#    if app.plot:
#        nk=int(round(sqrt(app.BZ.rank['k'])))
#        plt.title(engine.name.full+'_BC')
#        plt.axis('equal')
#        plt.colorbar(plt.pcolormesh(buff[:,0].reshape((nk,nk)),buff[:,1].reshape((nk,nk)),buff[:,2].reshape((nk,nk))))
#        if app.show:
#            plt.show()
#        else:
#            plt.savefig(engine.dout+'/'+engine.name.full+'_BC.png')
#        plt.close()
