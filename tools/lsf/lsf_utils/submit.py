import os
import shutil
import shlex
import subprocess

ENV_VARIABLES={
    "PrgEnv-intel":{"I_MPI_FABRICS_LIST":"tcp", 
                    "I_MPI_HYDRA_BOOTSTRAP":"lsf",
                    "I_MPI_HYDRA_BRANCH_COUNT":"-1",
                    "I_MPI_LSF_USE_COLLECTIVE_LAUNCH":"1"
                    }
}

'''
Base class for submitting jobs set up through a tcsh script. 
Base behaviour is to simply setup the environments in a shell script 
and launch that way.
'''
class TcshJob(object):
    def __init__(self):
        self._script_lines=[] #list the manage script lines to make it easier for different submission modes to merge 
        self.run_path = os.environ["PWD"] #Run the job in this directory
        self.stdout="tcsh_job.stdout"
        self.stderr="tcsh_job.stderr"
        self.cores=1
        self.command=None
        self.using_ncsu_mpi=False
        #environment setup
        self.load_modules = None
        self.using_conda_env=None
        self._generate_name="tcsh_job_script.csh"
        if 'PYTHON_ENV' in os.environ:
            self.python_env = os.environ['PYTHON_ENV']
        else:
            self.python_env = ''
            
    def set_env(self,l,m):#handle setting environment variables for certain modules
        if m not in ENV_VARIABLES: return
        for variable in ENV_VARIABLES[m]:
            l.append('setenv {} {}'.format(variable,ENV_VARIABLES[m][variable]))

    def final_execution_command(self):
        #return a string that will serve as the final executing command
        if self.using_ncsu_mpi:
            return "mpirun -n {} {}".format(self.cores,self.command)
        else:
            return "{}".format(self.command)
    
    def generate(self):
        #generate the main body of the script to execute the job
        self._script_lines.append("#!/bin/tcsh")
        #self._script_lines.append("cd {}".format(self.run_path))
        if len(self.load_modules) > 0:
            self._script_lines.append('module load {}'.format(" ".join(self.load_modules)))
            for m in self.load_modules: self.set_env(self._script_lines,m) #set additional env variables for certain modules
        if len(self.python_env) > 0:
            self._script_lines.append('source {}'.format(self.python_env))                
        if self.using_conda_env:
            self._script_lines.append('conda activate {}'.format(self.using_conda_env))               
        self._script_lines.append(self.final_execution_command())

    def submission_command(self,submission_script_name):
        return "(tcsh {} > {})>& {}".format(submission_script_name,self.stdout,self.stderr)

    def dump_script_lines(self):
        #dump script lines to the job script we are generating
        with open(self._generate_name,"w+") as f:
            for l in self._script_lines:
                f.write("{}\n".format(l))
                
    def submit(self, no=False):
        self._script_lines=[] #clear out script lines before writing new ones
        assert self.command!=None
        current_path = os.environ["PWD"]
        self.generate()
        self.dump_script_lines()
        shutil.copy(self._generate_name,self.run_path)
        os.chdir(self.run_path)
        if no == False:
            c=self.submission_command(self._generate_name)
            subprocess.Popen(c,shell=True,stdin=None,stdout=None,stderr=None)
        os.chdir(current_path)

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
    def using_ncsu_mpi(self):
        return self._using_ncsu_mpi
    @using_ncsu_mpi.setter
    def using_ncsu_mpi(self,x):
        self._using_ncsu_mpi=x

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
    
    @property 
    def using_conda_env(self):
        return self._using_conda_env
    @using_conda_env.setter
    def using_conda_env(self,env_path):
        self._using_conda_env=env_path

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
    def cores(self):
        return self._cores
    @cores.setter
    def cores(self,n):
        self._cores=n

'''
Class that helps with submitting jobs to an LSF scheduler, can set mostly any option you want
once ready call submit() and the job will be run with the corresponding parameters
at the given run path. It is up to the specific experiment to set these properties.
'''
class LSFJob(TcshJob):
    def __init__(self):
        TcshJob.__init__(self)
        self.memory=None  #8 GB memory
        self.time = 10 #10 Minutes
        self.queue="tuck"
        self.avoid_hosts=[]
        self.stdout="lsfjob.stdout"
        self.stderr="lsfjob.stderr"
        self.job_name="job"
        self.exclusive=False
        self._generate_name="hpc_lsfjob.csh"

    
    def submission_command(self,submission_script_name):
        return "bsub < {}".format(submission_script_name)


    def final_execution_command(self):
        #return a string that will serve as the final executing command, allows for command line specialization based on execution backend
        if self.using_ncsu_mpi:
            return  "mpirun {}".format(self.command)
        else:
            return "{}".format(self.command)
    
    def generate(self):
        TcshJob.generate(self) # generate the tcsh portion of the script
        lsf_script_lines = [] 
        lsf_script_lines.append("#BSUB -n " + str(self.cores))
        lsf_script_lines.append("#BSUB -W {}:00".format(self.time))
        if self.exclusive: lsf_script_lines.append("#BSUB -x")
        if self.one_host: lsf_script_lines.append("#BSUB -R \"span[hosts=1]\"") #this is to make sure all processes/thread sit at one host
        if self.memory: lsf_script_lines.append("#BSUB -R \"rusage[mem={}GB]\"".format(self.memory))
        if len(self.avoid_hosts)>0:
            #avoid problematic hosts
            for h in self.avoid_hosts:
                lsf_script_lines.append("#BSUB -R \"hname != {} \"".format(h))
        lsf_script_lines.append("#BSUB -q "+self.queue)
        lsf_script_lines.append("#BSUB -J " + self.job_name)
        lsf_script_lines.append("#BSUB -o {}.%J".format(self.stdout))
        lsf_script_lines.append("#BSUB -e {}.%J".format(self.stderr))
        self._script_lines=[self._script_lines[0]] + lsf_script_lines + self._script_lines[1::] #insert bsub header 
        
    @property
    def job(self):
        return self._job
    @job.setter
    def job(self,j):
        self._job=j
        
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
    def exclusive(self):
        return self._exclusive
    @exclusive.setter
    def exclusive(self,e):
        self._exclusive=e

    @property
    def one_host(self):
        return self._one_host
    @one_host.setter
    def one_host(self,e):
        self._one_host=e

    @property
    def avoid_hosts(self):
        return self._avoid_hosts
    @avoid_hosts.setter
    def avoid_hosts(self,hosts):
        if not type(hosts) is list:
            h=[hosts]
        else:
            h=hosts
        self._avoid_hosts=h

