import os
import shutil

'''
Class that helps with submitting jobs in python, can set mostly any option you want
once ready call submit() and the job will be run with the corresponding parameters
at the given run path. It is up to the specific experiment to set these properties.
'''
class LSFJob(object):
    def __init__(self):
        self.cores=1 #1 core
        self.memory=8  #8 GB memory
        self.time = 10 #10 Minutes
        self.run_path = os.environ["PWD"] #Run the job in this directory
        self.queue="tuck"
        self.command=None
        self.stdout="lsfjob.stdout"
        self.stderr="lsfjob.stderr"
        self.load_modules = ['conda']
        if 'PYTHON_ENV' in os.environ:
            self.python_env = os.environ['PYTHON_ENV']
        else:
            self.python_env = ''
        self.job="job"
        self.exclusive=False

    def generate(self):
        #generate the actual file
        generate_name = "hpc_lsfjob.csh"
        with open(generate_name,"w+") as f:
            f.write("#BSUB -n " + str(self.cores))
            f.write("\n#BSUB -W {}:00".format(self.time))
            if self.exclusive: f.write("\n#BSUB -x")
            #f.write("\n#BSUB -R span[hosts=1]")
            f.write("\n#BSUB -R \"rusage[mem={}GB]\"".format(self.memory))
            #f.write("\n#BSUB -M {}GB!".format(str(self.memory)))
            f.write("\n#BSUB -q "+self.queue)
            f.write("\n#BSUB -J " + self.job)
            f.write("\n#BSUB -o {}.%J".format(self.stdout))
            f.write("\n#BSUB -e {}.%J".format(self.stderr))
            if len(self.load_modules) > 0:
                f.write('\nmodule load {}'.format(" ".join(self.load_modules)))
            if len(self.python_env) > 0:
                f.write('\nsource {}'.format(self.python_env))
            f.write("\n {}".format(self.command))
        return generate_name
    def submit(self, no=False):
        assert self.command!=None
        current_path = os.environ["PWD"]
        generate_name = self.generate()
        shutil.copy(generate_name,self.run_path)
        os.chdir(self.run_path)
        if no == False:
            os.system('bsub <' +generate_name)
        os.chdir(current_path)

    @property
    def job(self):
        return self._job
    @job.setter
    def job(self,j):
        self._job=j
        
    @property
    def cores(self):
        return self._cores
    @cores.setter
    def cores(self,n):
        self._cores=n

    @property
    def memory(self):
        return self._memory
    @memory.setter
    def memory(self,n):
        self._memory=n

        
    @property
    def time(self):
        return self._time
    @time.setter
    def time(self,t):
        self._time = t

    @property
    def queue(self):
        return self._queue
    @queue.setter
    def queue(self,q:str):
        self._queue=q

    @property
    def stdout(self):
        return self._stdout
    @stdout.setter
    def stdout(self,s:str):
        self._stdout=s

    @property
    def stderr(self):
        return self._stderr
    @stderr.setter
    def stderr(self,s:str):
        self._stderr=s

    @property
    def command(self):
        return self._executeable_path
    @command.setter
    def command(self,s:str):
        self._executeable_path=s

    @property
    def run_path(self):
        return self._run_path
    @run_path.setter
    def run_path(self,s:str):
        self._run_path=s
       
    @property
    def exclusive(self):
        return self._exclusive
    @exclusive.setter
    def exclusive(self,e):
        self._exclusive=e

    @property
    def load_modules(self):
        return self._load_modules
    @load_modules.setter
    def load_modules(self,l):
        self._load_modules=l

    @property
    def python_env(self):
        return self._python_env
    @python_env.setter
    def python_env(self,p):
        assert type(p)==str and "python_env should be a string."
        #assert os.path.exists(p) and "python_env should exist."
        self._python_env=p
