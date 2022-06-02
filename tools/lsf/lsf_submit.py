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
        self.job="job"
        self.exclusive=False

    def generate(self):
        #generate the actual file
        generate_name = "hpc_lsfjob.csh"
        with open(generate_name,"w+") as f:
            f.write("#BSUB -n " + str(self.cores))
            f.write("\n#BSUB -W " +str(self.time))
            if self.exclusive: f.write("\n#BSUB -x")
            f.write("\n#BSUB -R span[hosts=1]")
            f.write("\n#BSUB -R \"rusage[mem={}GB]\"".format(self.memory))
            f.write("\n#BSUB -M {}GB!".format(str(self.memory)))
            f.write("\n#BSUB -q "+self.queue)
            f.write("\n#BSUB -J " + self.job)
            f.write("\n#BSUB -o {}.%J".format(self.stdout))
            f.write("\n#BSUB -e {}.%J".format(self.stderr))
            f.write("\n {}".format(self.command))
        return generate_name
    def submit(self):
        assert self.command!=None
        current_path = os.environ["PWD"]
        generate_name = self.generate()
        shutil.copy(generate_name,self.run_path)
        os.chdir(self.run_path)
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
