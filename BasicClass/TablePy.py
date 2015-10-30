from IndexPy import *
class Table(dict):
    '''
    The Table class provides the methods to get an index from its sequence number or vice versa.
    '''
    def __init__(self,indices=[],dict={},f=None):
        '''
        If the function f is assigned, the sequence number corresponding to a index will be set by f. Otherwise it equals the index's order in the indices list. An ordinary dictionary can also be used to initialize an instance.
        '''
        for i,v in enumerate(indices):
            if f==None:
                self[v]=i
            else:
                self[v]=f(v)
        self.update(dict)
           
def union(**kwargs):
    '''
    Define the union of index-sequence tables.
    '''
    result=Table()
    sum=0
    for k,v in kwargs.iteritems():
        if isinstance(v,Table):
            result[k]=v
            count=0
            for k1,v1 in v.iteritems():
                if not isinstance(v1,Table): 
                    result[k1]=v1+sum
                    count+=1
            sum+=count
    return result

def subset(table,mask):
    '''
    Define a certain subset of a sequence table according to the mask function.
    '''
    result=Table()
    result['_superset']=table
    for k,v in table.iteritems():
        if not isinstance(v,Table) and mask(k):
            result[k]=v
    buff={}
    for i,k in enumerate(sorted([key for key in result.keys() if not isinstance(result[key],Table)],key=result.get)):
        buff[k]=i
    result.update(buff)
    return result
    
def reverse_table(table):
    '''
    Generate the sequence-index table for a reversed lookup.
    '''
    result=Table()
    for k,v in table.iteritems():
        if isinstance(v,Table):
            result[k]=reverse_table(v)
        else:
            result[v]=k
    return result
