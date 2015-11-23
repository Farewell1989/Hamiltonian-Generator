'''
Table.
'''
class Table(dict):
    '''
    This class provides the methods to get an index from its sequence number or vice versa.
    '''
    def __init__(self,indices=[],dict={},f=None):
        '''
        Constructor.
        Parameters:
            indices: list of Index
                The indices that need to be mapped to sequences.
            dict: dict, optional
                An already constructed index-sequence table.
            f: function, optional
                The function used to map an index to a sequence.
                If it is None, the order of the index in indices will be used as its sequence number.
        '''
        for i,v in enumerate(indices):
            if f==None:
                self[v]=i
            else:
                self[v]=f(v)
        self.update(dict)
           
def union(**kwargs):
    '''
    This function returns the union of index-sequence tables.
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
    This function returns a certain subset of an index-sequence table according to the mask function.
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
    This function returns the sequence-index table for a reversed lookup.
    '''
    result=Table()
    for k,v in table.iteritems():
        if isinstance(v,Table):
            result[k]=reverse_table(v)
        else:
            result[v]=k
    return result
