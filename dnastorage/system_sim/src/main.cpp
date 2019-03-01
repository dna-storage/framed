#include "system.h"
#include "decoder.h"
#include "prep.h"
#include "generator.h"
#include "scheduler.h"
#include "storage_system.h"
#include "sequencer.h"
#include "utlist.h"
#include "parameters.h"
#include "stats.h"
#include "json/include/rapidjson/document.h"
#include "json/include/rapidjson/filereadstream.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

using namespace rapidjson;

int main(int argc, char** argv){
  system_sim_t* system_simulator;
  system_sim_params_t sim_params;
  FILE* config;
  std::string phase_log;
  std::string trace_file;
  std::string stats_log;
  
  if(argc<2){
    printf("Provide a configuration file\n");
    return 1;
  }

  config=fopen(argv[1],"rb");
  if(config==NULL){
    printf("Cannot open config file\n");
    return 1;
  }

  char readBuffer[65536];
  FileReadStream s(config,readBuffer,sizeof(readBuffer));

  Document config_doc;
  config_doc.ParseStream(s);

  //go through the parsed json and get parameters
  
  



  
  
  return 0;

}
