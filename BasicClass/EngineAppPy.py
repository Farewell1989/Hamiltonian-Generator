import time
class Engine(object):
    '''
    Class Engine is the parent class for all classes of Hamiltonian-based algorithms. It has the following attributes:
    1) din: the directory where its needed data comes from;
    2) dout: the directory where its generated data writes to;
    3) waiting_list: the names of apps waiting to be run;
    4) apps: a dict containing all apps associated with this engine. 
    Note that not all apps are in the waiting list because some of them are actually called by others and do not need automatic running.
    '''

    def __new__(cls,*arg,**karg):
        result=object.__new__(cls,*arg,**karg)
        if 'din' in karg:
            result.din=karg['din']
        else:
            result.din='.'
        if 'dout' in karg:
            result.dout=karg['dout']
        else:
            result.dout='.'
        result.waiting_list=[]
        result.apps={}
        return result

    def __init__(self,*arg,**karg):
        pass

    def addapps(self,name=None,app=None):
        '''
        Add an app to the app dict with its key equals the parameter name if assigned or otherwise the app's type name. Unless the parameter name is assigned, it will not go into the waiting list.
        '''
        if name!=None:
            self.apps[name]=app
            self.waiting_list.append(name)
        else:
            self.apps[app.__class__.__name__]=app

    def runapps(self,name=None,clock=False):
        '''
        Run a specific app if the parameter name if assigned otherwise run all the apps whose names are contained in the waiting list in order. In the former case, an extra parameter clock can be used to determine whether the engine records the time the app consumes. In the latter case, the consumed time is recorded for every app in the waiting list respectively.
        '''
        if name!=None:
            if clock :
                stime=time.time()
                self.apps[name].run(self,self.apps[name])
                etime=time.time()
                print 'App '+name+': time consumed '+str(etime-stime)+'s.'
            else:
                self.apps[name].run(self,self.apps[name])
        else:
            for name in self.waiting_list:
                stime=time.time()
                self.apps[name].run(self,self.apps[name])
                etime=time.time()
                print 'App '+name+': time consumed '+str(etime-stime)+'s.'

class App(object):
    '''
    Class App is the parent class for all classes of specific tasks based on the Engine. It has the following attributes:
    1) plot: whether the results are to be plotted;
    2) show: whether the plotted graph is to be shown;
    3) parallel: whether the calculating process is parallel;
    4) np: the number of processes used in parallel computing, and 0 means the available maximum;
    5) save_data: whether the generated data is to be saved on the hard drive;
    6) run: the function called by the engine to carry out the tasks, which should be complemented by the inherited engine.
    '''
    
    def __new__(cls,*arg,**karg):
        result=object.__new__(cls,*arg,**karg)
        if 'plot' in karg:
            result.plot=karg['plot']
        else:
            result.plot=True
        if 'show' in karg:
            result.show=karg['show']
        else:
            result.show=True
        if 'parallel' in karg:
            result.parallel=karg['parallel']
        else:
            result.parallel=False
        if 'np' in karg:
            result.np=karg['np']
        else:
            result.np=0
        if 'save_data' in karg:
            result.save_data=karg['save_data']
        else:
            result.save_data=True
        if 'run' in karg and callable(karg['run']):
            result.run=karg['run']
        return result

    def __init__(self,*arg,**karg):
        pass
