#include "dna_storage_attributes.h"
#include "decoder.h"
#include "prep.h"
#include "generator.h"
#include "scheduler.h"
#include "storage_system.h"
#include "sequencer.h"
#include "utlist.h"

#include<stdlib.h>
#include<stdio.h>


//function to strip the configuration file




int main(int argc, char** argv){
  system_sim_t* system_simulator;
  FILE* config_file;
  FILE* trace_file;
  //variables to be sent to the 

  config_file=fopen(argv[1],"r");
  if(config_file==NULL){
    perror("Missing config file");
    return 1;
  }
  trace_file=fopen(argv[1],"r");
  if(trace_file==NULL){
    perror("Missing trace file");
    return 1;
  }

  //instantiate the system object
  //system_simulator=new system();
  //parse the configuration file

  

  //simulate the system
  //system_simulator->simulate();

  //done with simulation
  //FIX ME: Need to print out statistics
  //delete system_simulator;
  return 0;

}
