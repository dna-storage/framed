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
  FILE* phase;
  FILE* stats;
  FILE* trace;
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
  sim_params.num_preps=config_doc["prep"]["number"].GetUint64();
  sim_params.prep_channels=config_doc["prep"]["channels"].GetUint64();
  sim_params.prep_time=config_doc["prep"]["time"].GetUint64();
  sim_params.num_sequencers=config_doc["sequencer"]["number"].GetUint64();
  sim_params.max_strands_sequencer=config_doc["sequencer"]["max_strands"].GetUint64();
  sim_params.seq_channels=config_doc["sequencer"]["channels"].GetUint64();
  sim_params.seq_time=config_doc["sequencer"]["time"].GetUint64();
  sim_params.seq_efficiency=config_doc["sequencer"]["efficiency"].GetDouble();
  sim_params.sequencer_timeout=config_doc["sequencer"]["timeout"].GetUint64();
  sim_params.num_decoders=config_doc["decoder"]["number"].GetUint64();
  sim_params.dec_time=config_doc["decoder"]["time"].GetUint64();
  sim_params.average_pool_capacity=config_doc["storage"]["avg_pool_capacity"].GetUint64();
  sim_params.number_pools=config_doc["storage"]["number_pools"].GetUint64();
  sim_params.bytes_per_strand=config_doc["storage"]["bytes_per_strand"].GetDouble();
  sim_params.pool_write_time=config_doc["storage"]["pool_write_time"].GetUint64();
  sim_params.pool_wait_time=config_doc["storage"]["pool_wait_time"].GetUint64();
  sim_params.pool_copies=config_doc["storage"]["pool_copies"].GetUint64();
  sim_params.number_reads=config_doc["storage"]["number_reads"].GetUint64();
  sim_params.max_file_size=config_doc["storage"]["max_file_size"].GetUint64();
  sim_params.min_file_size=config_doc["storage"]["min_file_size"].GetUint64();
  sim_params.sequencing_depth=config_doc["storage"]["sequencing_depth"].GetUint64();
  sim_params.seed=config_doc["generator"]["seed"].GetInt();
  sim_params.rate=config_doc["generator"]["rate"].GetDouble();
  sim_params.batch_size=config_doc["scheduler"]["batch_size"].GetUint64();
  sim_params.window_size=config_doc["system"]["window_size"].GetUint64();
  sim_params.prep_seq_buffer_size=config_doc["system"]["prep_seq_buffer_size"].GetUint64();
  sim_params.seq_dec_buffer_size=config_doc["system"]["seq_dec_buffer_size"].GetUint64();
  sim_params.sim_time=config_doc["system"]["simulation_time"].GetUint64();
  phase_log=config_doc["phase_log"].GetString();
  if(!config_doc["trace_file"].IsNull())trace_file=config_doc["trace_file"].GetString();
  stats_log=config_doc["stats_log"].GetString();

  if(!config_doc["trace_file"].IsNull()) trace=fopen(trace_file.c_str(),"r");
  else trace=NULL;

  phase=fopen(phase_log.c_str(),"w");
  if(phase==NULL) {
    printf("phase log open error\n");
    exit(1);
  }
  stats=fopen(stats_log.c_str(),"w");
  if(stats==NULL){
    printf("stats log open error\n");
    exit(1);
  }
  
  sim_params.trace_file=trace;
  sim_params.phase_log=phase;
  sim_params.stats_log=stats;

  system_simulator=new system_sim_t(sim_params);

  system_simulator->simulate(); //run the simulator
  
  delete system_simulator; //delete the simulator to free memory

    
  return 0;

}
