FROM centos:7

RUN groupadd -g 699 appuser && \
    useradd -m -r -u 699 -g appuser appuser

USER root

WORKDIR /home/appuser/framed

LABEL maintainer="kevinvolkel95@gmail.com"


SHELL [ "/bin/bash", "--login", "-c" ]

RUN yum update -y \
 && yum install git zip wget bash tcsh python3-pip python3 gcc-c++ which -y \
 && yum groupinstall "Development Tools" -y


RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.2-linux-x86_64.tar.gz \
 && tar xvf julia-1.6.2-linux-x86_64.tar.gz \
 && mkdir /opt/julia \
 && mv julia-1.6.2 /opt/julia \
 && rm julia-1.6.2-linux-x86_64.tar.gz


RUN  wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.5.tar.gz \
 && tar xvf openmpi-4.1.5.tar.gz \
 && cd openmpi-4.1.5 \
 && ./configure --prefix /usr/lib64/openmpi-4.1.5 \
 && make all \
 && make install \
 && make -j 4 


ENV CONDA_DIR="/home/appuser/conda"
ENV PATH="$PATH:/opt/julia/julia-1.6.2/bin:$CONDA_DIR/bin:/usr/lib64/openmpi-4.1.5/bin"
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/lib64/openmpi-4.1.5/lib"


RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh \
 && /bin/bash ~/miniconda.sh -b -p /home/appuser/conda \ 
 && conda init bash


#still need to install FrameD itself, along with schwimbadd submodule

RUN chown appuser /home/appuser/framed

USER appuser

RUN touch ~/.tcshrc \
  && touch ~/.bashrc \
  && conda init --user tcsh \
  && conda init --user bash 


COPY ./dnastorage.yml ./dnastorage.yml
COPY ./requirements.txt ./requirements.txt

RUN conda env create -f dnastorage.yml


ENV CC="/usr/lib64/openmpi-4.1.5/bin/mpicc"
ENV MPICC="/usr/lib64/openmpi-4.1.5/bin/mpicc"

RUN conda activate dnastorage \
&& pip install -r requirements.txt


RUN conda activate dnastorage \
  && python -c "import julia; julia.install()" 


COPY --chown=appuser . . 

RUN conda activate dnastorage \
&& tcsh ./init.csh -no-env \
&& make develop \
&& cd schwimmbad \
&& pip install . \
&& cd ../ \
&& conda config --set auto_activate_base false

  
ENV PATH="$PATH:/usr/lib64/openmpi/bin"

ENV HOME="/home/appuser"
ENV DNASTORAGE_HOME="$HOME/framed"
ENV DNASTORAGE_TOOLS="$HOME/framed/tools"
ENV DNASTORAGE_LSF="$HOME/framed/tools/lsf"
ENV DNArSimPath="$HOME/framed/dnastorage/fi/DNArSim"
ENV FRAMED_CONFIGS="$HOME/framed/examples"
ENV FRAMED_CONDA="dnastorage"
ENV DNArSimData="$HOME/framed/probEdit"
ENV FRAMED_IMAGE_FILES="$HOME/framed/test_files"
ENV OMPI_MCA_btl_vader_single_copy_mechanism=none


ENTRYPOINT ["conda","run","--no-capture-output","-n","dnastorage","/bin/bash"]






 


 
