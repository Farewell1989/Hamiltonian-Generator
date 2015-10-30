from OperatorPy import *
from BasisEPy import *
from scipy.sparse import *
from numba import jit
def opt_rep(operator,basis,transpose=False,dtype=complex128):
    '''
    The opt_rep function returns the sparse matrix representation of an operator on the occupation number basis. Three kinds of operators are supported:
    1) The single particle operators. Since a single creation or annihilation operator changes the number of particles, the quantum state usually lives in a different Hilbert space after the operator's operation on it. Therefore in this case, the parameter basis is a list of two instances of BasisE. 
    2) The quadratic operators. Any kind of quadratic terms are supported, e.g., the hopping terms, onsite terms and pairing terms. In this case, only one basis is needed. Note that for pairing terms, only those basis with mode="EG"
    should be used because of the non-conservation of particle numbers. 
    3) The quartic operators. For simplicity, only the normal ordered quartic operators with two creation terms and two annihilation terms are supported.
    All of those operators' representations are generated in the real space. If the parameter transpose is set to be true, the returned sparse matrix is stored in the csr form, otherwise in the csc form. The difference is just a transpose.
    '''
    if operator.rank==1:
        if len(basis)==2:
            return opt_rep_1(operator,basis[0],basis[1],transpose,dtype=dtype)
        else:
            raise ValueError("Opt_rep error: when the operator's rank is 1 there must be two groups of basis.")
    elif operator.rank==2:
        return opt_rep_2(operator,basis,transpose,dtype=dtype)
    elif operator.rank==4:
        if operator.is_normal_ordered():
            return opt_rep_4(operator,basis,transpose,dtype=dtype)
        else:
            raise ValueError("Opt_rep error: when the operator's rank is 4, it must be normal ordered.")
    else:
        raise ValueError("Opt_rep error: only operators with rank=1,2,4 are supported.")

def opt_rep_1(operator,basis1,basis2,transpose,dtype=complex128):
    nbasis1=basis1.nbasis
    nbasis2=basis2.nbasis
    data=zeros(nbasis1,dtype=dtype)
    indices=zeros(nbasis1,dtype=int32)
    indptr=zeros(nbasis1+1,dtype=int32)
    seq=operator.seqs[0]
    basis_table1=basis1.basis_table
    basis_table2=basis2.basis_table
    if operator.indices[0].nambu==CREATION:
        opt_rep_1_1(data,indices,indptr,nbasis1,nbasis2,basis_table1,basis_table2,seq)
    else:
        opt_rep_1_0(data,indices,indptr,nbasis1,nbasis2,basis_table1,basis_table2,seq)
    if transpose==False:
        return csr_matrix((data*operator.value,indices,indptr),shape=(nbasis1,nbasis2))
    else:
        return csr_matrix((data*operator.value,indices,indptr),shape=(nbasis1,nbasis2)).T

@jit
def opt_rep_1_1(data,indices,indptr,nbasis1,nbasis2,basis_table1,basis_table2,seq):
    ndata=0
    for i in xrange(nbasis1):
        indptr[i]=ndata
        rep=basis_rep(i,basis_table1)
        if not rep&1<<seq:
            nsign=0
            for j in xrange(seq):
                if rep&1<<j: nsign+=1
            rep=rep|1<<seq
            indices[ndata]=seq_basis(rep,basis_table2)
            data[ndata]=(-1)**nsign
            ndata+=1
    indptr[nbasis1]=ndata

@jit
def opt_rep_1_0(data,indices,indptr,nbasis1,nbasis2,basis_table1,basis_table2,seq):
    ndata=0
    for i in xrange(nbasis1):
        indptr[i]=ndata
        rep=basis_rep(i,basis_table1)
        if rep&1<<seq:
            nsign=0
            for j in xrange(seq):
                if rep&1<<j: nsign+=1
            rep=rep&~(1<<seq)
            indices[ndata]=seq_basis(rep,basis_table2)
            data[ndata]=(-1)**nsign
            ndata+=1
    indptr[nbasis1]=ndata

def opt_rep_2(operator,basis,transpose,dtype=complex128):
    nbasis=basis.nbasis
    data=zeros(nbasis,dtype=dtype)
    indices=zeros(nbasis,dtype=int32)
    indptr=zeros(nbasis+1,dtype=int32)
    seq1=operator.seqs[0]
    seq2=operator.seqs[1]
    basis_table=basis.basis_table
    if operator.indices[0].nambu==CREATION and operator.indices[1].nambu==ANNIHILATION:
        opt_rep_2_10(data,indices,indptr,nbasis,basis_table,seq1,seq2)
    elif operator.indices[0].nambu==ANNIHILATION and operator.indices[1].nambu==CREATION:
        opt_rep_2_01(data,indices,indptr,nbasis,basis_table,seq1,seq2)
    elif opetator.indices[0].nambu==ANNIHILATION and operator.indices[1].nambu==ANNIHILATION:
        opt_rep_2_00(data,indices,indptr,nbasis,basis_table,seq1,seq2)
    else:
        opt_rep_2_11(data,indices,indptr,nbasis,basis_table,seq1,seq2)
    if transpose==False:
        return csr_matrix((data*operator.value,indices,indptr),shape=(nbasis,nbasis))
    else:
        return csr_matrix((data*operator.value,indices,indptr),shape=(nbasis,nbasis)).T

@jit
def opt_rep_2_10(data,indices,indptr,nbasis,basis_table,seq1,seq2):
    ndata=0
    for i in xrange(nbasis):
        indptr[i]=ndata
        rep=basis_rep(i,basis_table)
        if rep&1<<seq2:
            rep1=rep&~(1<<seq2)
            if not rep1&(1<<seq1):
                nsign=0
                for j in xrange(seq2):
                    if rep&1<<j: nsign+=1
                for j in xrange(seq1):
                    if rep1&1<<j: nsign+=1
                rep2=rep1|1<<seq1
                indices[ndata]=seq_basis(rep2,basis_table)
                data[ndata]=(-1)**nsign
                ndata+=1
    indptr[nbasis]=ndata

@jit
def opt_rep_2_01(data,indices,indptr,nbasis,basis_table,seq1,seq2):
    ndata=0
    for i in xrange(nbasis):
        indptr[i]=ndata
        rep=basis_rep(i,basis_table)
        if not rep&1<<seq2:
            rep1=rep|1<<seq2
            if rep1&1<<seq1:
                nsign=0
                for j in xrange(seq2):
                    if rep&1<<j: nsign+=1
                for j in xrange(seq1):
                    if rep1&1<<j: nsign+=1
                rep2=rep1&~(1<<seq1)
                indices[ndata]=seq_basis(rep2,basis_table)
                data[ndata]=(-1)**nsign
                ndata+=1
    indptr[nbasis]=ndata

@jit
def opt_rep_2_00(data,indices,indptr,nbasis,basis_table,seq1,seq2):
    ndata=0
    for i in xrange(nbasis):
        indptr[i]=ndata
        rep=basis_rep(i,basis_table)
        if rep&1<<seq2:
            rep1=rep&~(1<<seq2)
            if rep1&1<<seq1:
                nsign=0
                for j in xrange(seq2):
                    if rep&1<<j: nsign+=1
                for j in xrange(seq1):
                    if rep1&1<<j: nsign+=1
                rep2=rep1&~(1<<seq1)
                indices[ndata]=seq_basis(rep2,basis_table)
                data[ndata]=(-1)**nsign
                ndata+=1
    indptr[nbasis]=ndata

@jit
def opt_rep_2_11(data,indices,indptr,nbasis,basis_table,seq1,seq2):
    ndata=0
    for i in xrange(nbasis):
        indptr[i]=ndata
        rep=basis_rep(i,basis_table)
        if not rep&1<<seq2:
            rep1=rep|1<<seq2
            if not rep1&1<<seq1:
                nsign=0
                for j in xrange(seq2):
                    if rep&1<<j: nsign+=1
                for j in xrange(seq1):
                    if rep1&1<<j: nsign+=1
                rep2=rep1|1<<seq1
                indices[ndata]=seq_basis(rep2,basis_table)
                data[ndata]=(-1)**nsign
                ndata+=1
    indptr[nbasis]=ndata

def opt_rep_4(operator,basis,transpose,dtype=complex128):
    nbasis=basis.nbasis
    data=zeros(nbasis,dtype=dtype)
    indices=zeros(nbasis,dtype=int32)
    indptr=zeros(nbasis+1,dtype=int32)
    seq1=operator.seqs[0]
    seq2=operator.seqs[1]
    seq3=operator.seqs[2]
    seq4=operator.seqs[3]
    basis_table=basis.basis_table
    opt_rep_4_1100(data,indices,indptr,nbasis,basis_table,seq1,seq2,seq3,seq4)
    if transpose==False:
        return csr_matrix((data*operator.value,indices,indptr),shape=(nbasis,nbasis))
    else:
        return csr_matrix((data*operator.value,indices,indptr),shape=(nbasis,nbasis)).T

@jit
def opt_rep_4_1100(data,indices,indptr,nbasis,basis_table,seq1,seq2,seq3,seq4):
    ndata=0
    for i in xrange(nbasis):
        indptr[i]=ndata
        rep=basis_rep(i,basis_table)
        if rep&1<<seq4:
            rep1=rep&~(1<<seq4)
            if rep1&1<<seq3:
                rep2=rep1&~(1<<seq3)
                if not rep2&1<<seq2:
                    rep3=rep2|1<<seq2
                    if not rep3&1<<seq1:
                        nsign=0
                        for j in xrange(seq4):
                            if rep&1<<j: nsign+=1
                        for j in xrange(seq3):
                            if rep1&1<<j: nsign+=1
                        for j in xrange(seq2):
                            if rep2&1<<j: nsign+=1
                        for j in xrange(seq1):
                            if rep3&1<<j: nsign+=1
                        rep4=rep3|1<<seq1
                        indices[ndata]=seq_basis(rep4,basis_table)
                        data[ndata]=(-1)**nsign
                        ndata+=1
    indptr[nbasis]=ndata
