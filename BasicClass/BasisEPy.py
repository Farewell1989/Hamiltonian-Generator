from numpy import *
from math import factorial
from itertools import combinations
from numba import jit

class BasisE:
    '''
    Class BasisE deals with the occupation-number basis of electron systems and provides a unified description of the three often-encountered cases, i.e. particle-non-conserved systems, particle-conserved systems and spin-conserved systems. Its attributes are as follows:
    1) basis_type: a string to distinguish the type of the three kinds of fore mentioned systems;
    2) nstate: an array containing the numbers of states;
    3) nparticle: an array containing the numbers of particles;
    4) basis_table: the table of allowed binary basis of the Hilbert space;
    5) nbasis: the dimension of the Hilbert space.
    '''
    
    def __init__(self,tuple=(),up=(),down=(),nstate=0,dtype=int64):
        '''
        There are three ways to initialize a BasisE instance:
        1) Assign the parameter nstate only, which will generate the basis for a particle-non-conserved system with nstate available generalized orbitals(the site and spin index are interpreted as generalized orbital index). The attribute basis_type will be "EG".
        2) Assign the parameter tuple only, which will generate the basis for a particle-conserved system with tuple[0] generalized orbitals occupied by tuple[1] electrons. The attribute basis_type will be "EP".
        3) Assign the up and down tuples both, which will generate the basis for a spin-conserved system with up[0] orbitals occupied by up[1] spin-up electrons and down[0] orbitals occupied by down[1] spin-down electrons. The attribute basis_type will be "ES".
        Note: if more than needed parameters to generate a certain kind a basis are assigned, the __init__ method obeys the following priority to create the instance: "EP" > "ES" > "EG".
        '''
        if len(tuple)==2:
            self.basis_type="EP"
            self.nstate=array(tuple[0])
            self.nparticle=array(tuple[1])
            self.basis_table=basis_table_ep(tuple[0],tuple[1],dtype=dtype)
            self.nbasis=len(self.basis_table)
        elif len(up)==2 and len(down)==2:
            self.basis_type="ES"
            self.nstate=array([up[0],down[0]])
            self.nparticle=array([up[1],down[1]])
            self.basis_table=basis_table_es(up,down,dtype=dtype)
            self.nbasis=len(self.basis_table)
        else:
            self.basis_type="EG"
            self.nstate=array(nstate)
            self.nparticle=array([])
            self.basis_table=array([])
            self.nbasis=2**nstate

    def __str__(self):
        '''
        Convert an instance to string.
        '''
        result=''
        if self.basis_type=='EG':
            for i in xrange(self.nbasis):
                result+=str(i)+': '+'{0:b}'.format(i)+'\n'
        else:
            for i,v in enumerate(self.basis_table):
                result+=str(i)+': '+'{0:b}'.format(v)+'\n'
        return result

def basis_table_ep(nstate,nparticle,dtype=int64):
    '''
    Generate the binary basis table with nstate orbitals occupied by nparticle electrons.
    '''
    result=zeros(factorial(nstate)/factorial(nparticle)/factorial(nstate-nparticle),dtype=dtype)
    buff=combinations(xrange(nstate),nparticle)
    for i,v in enumerate(buff):
        basis=0
        for num in v:
            basis+=(1<<num)
        result[i]=basis
    result.sort()
    return result

def basis_table_es(up,down,dtype=int64):
    '''
    Generate the binary basis table according to the up and down tuples.
    '''
    result=zeros(factorial(up[0])/factorial(up[1])/factorial(up[0]-up[1])*factorial(down[0])/factorial(down[1])/factorial(down[0]-down[1]),dtype=dtype)
    buff_up=list(combinations(xrange(up[0]),up[1]))
    buff_dn=list(combinations(xrange(down[0]),down[1]))
    count=0
    for vup in buff_up:
        buff=0
        for num in vup:
            buff+=(1<<num)
        buff=buff<<down[0]
        for vdn in buff_dn:
            basis=buff
            for num in vdn:
                basis+=(1<<num)
            result[count]=basis
            count+=1
    result.sort()
    return result

@jit
def basis_rep(seq_basis,basis_table):
    '''
    Find the binary basis representation whose sequence in basis_table is seq_basis.
    '''
    if len(basis_table)==0:
        return seq_basis
    else:
        return basis_table[seq_basis]

@jit
def seq_basis(basis_rep,basis_table):
    '''
    Find the basis sequence of basis_rep in basis_table.
    '''
    if len(basis_table)==0 :
        return basis_rep
    else:
        lb=0;ub=len(basis_table)
        result=(lb+ub)/2
        while basis_table[result]!=basis_rep:
            if basis_table[result]>basis_rep:
                ub=result
            else:
                lb=result
            if ub==lb: 
                error=True
                break
            result=(lb+ub)/2
        else:
            return result
        raise ValueError('Seq_basis error: the input basis_rep is not in the basis_table.')
