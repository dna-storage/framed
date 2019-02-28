/*
File: stats.cpp
Description: Member functions for the class that provides help with collecting statistics
Credit and Author: This code was adapted from Dr. Rotenberg's 721sim microarchitecture simulator.
*/

#include "stats.h"
#include "parameters.h"
#include <stdlib.h>
#include <cinttypes>
#include <cstring>
#include <map>
#include <cstdio>

stats_t::stats_t(FILE* stats_log, FILE* phase_log){
  this->set_log_files(stats_log,phase_log); //set the log files for the stat class

  /*
    Declare counters for information to be collected 

  */

  
  DECLARE_COUNTER(this, time_step, dna_storage); //counters to track timestep
  DECLARE_PHASE_COUNTER(this, time_step,dna_storage);

 
  DECLARE_COUNTER(this, finished_requests, dna_storage); //counters to track number of requests finished
  DECLARE_PHASE_COUNTER(this, finished_requests, dna_storage);

  DECLARE_COUNTER(this, data_decoded, dna_storage); //track the total amount of data decoded and provided back to users
  DECLARE_PHASE_COUNTER(this, data_decoded, dna_storage);


  DECLARE_COUNTER(this, pool_writes, dna_storage); //track the number of times a pool needs to be re-synthesized
  DECLARE_PHASE_COUNTER(this, pool_writes, dna_storage);

  DECLARE_COUNTER(this, total_latency, dnastorage); //track the total latency of all finished requests... used for tracking average latency
  DECLARE_PHASE_COUNTER(this, total_latency, dnastorage);

  
  
  /*
    Declare rates, rates should be a function of previously declared counters
  */
  
  DECLARE_RATE(this, IOPH, dna_storage, finished_requests, time_step, ((float)STEP_PER_HOUR)); //track the I/O operations per hour
  DECLARE_PHASE_RATE(this, IOPH, dna_storage, finished_requests, time_step, ((float)STEP_PER_HOUR));

  DECLARE_RATE(this, BW, dna_storage, data_decoded, time_step, ((float)FILE_UNIT*(float)STEP_PER_HOUR)); //track the raw data bandwidth: units -> Bytes/Hour
  DECLARE_PHASE_RATE(this, BW, dna_storage, data_decoded, time_step, ((float)FILE_UNIT*(float)STEP_PER_HOUR));

  
  DECLARE_RATE(this, writes_per_request, dna_storage, pool_writes, finished_requests,1.0); //track the rate of pool writes per finished request
  DECLARE_PHASE_RATE(this,writes_per_request, dna_storage, pool_writes, finished_requests,1.0);
  

  DECLARE_RATE(this, average_latency, dna_storage, total_latency, finished_requests,1.0); //track the average latency for a request in hours
  DECLARE_PHASE_RATE(this, average_latency, dna_storage, total_latency, finished_requests,1.0/((float)STEP_PER_HOUR));


  

  reset_counters();
  reset_phase_counters();
  set_phase_interval("time_step",100000);
  phase_id = 0;

}

void stats_t::set_log_files(FILE* _stats_log,FILE* _phase_log){
  this->stats_log = _stats_log;
  this->phase_log = _phase_log;
}

void stats_t::set_phase_interval(const char* name,uint64_t interval)
{
  std::strcpy(phase_counter_name,name);
  phase_interval = interval;
  fprintf(stderr,"Setting phase interval to %s = %lu\n",phase_counter_name,interval);
}

void stats_t::reset_counters(){
  std::map<std::string, counter_t*, ltstr>::iterator ctr_iter;
  for(ctr_iter = counter_map.begin();ctr_iter != counter_map.end(); ctr_iter++){
    ctr_iter->second->count = 0;
  }
}

void stats_t::reset_phase_counters(){
  std::map<std::string, counter_t*, ltstr>::iterator ctr_iter;
  for(ctr_iter = counter_map.begin();ctr_iter != counter_map.end(); ctr_iter++){
    ctr_iter->second->phase_count = 0;
  }
}

void stats_t::register_counter(const char* name, const char* hierarchy){

  counter_t* c    = new counter_t;
  c->count        = 0;
  c->phase_count  = 0;
  c->name         = new char[strlen(name)+1];
  c->hierarchy    = new char[strlen(hierarchy)+1];
  c->valid_phase_counter    = false;
  strcpy(c->name,name);
  strcpy(c->hierarchy,hierarchy);
  counter_map[name] = c;
  fprintf(stderr,"Counter name %s %s\n",name,hierarchy);
  fflush(0);
}

void stats_t::register_phase_counter(const char* name, const char* hierarchy){
  // If the counter has been declared, mark it as a phase counter
  if(counter_map.find(name) != counter_map.end()){
    counter_map[name]->valid_phase_counter = true;
  } 
  // If it does not exist, declare it and mark it as a phase counter
  else {
    counter_t* c    = new counter_t;
    c->name         = new char[strlen(name)+1];
    c->hierarchy    = new char[strlen(hierarchy)+1];
    c->valid_phase_counter  = true;
    strcpy(c->name,name);
    strcpy(c->hierarchy,hierarchy);
    counter_map[name] = c;
  }
}

void stats_t::register_rate(const char* name, const char* hierarchy, const char* numerator, const char* denominator, double multiplier){
  rate_t* r       = new rate_t;
  r->name         = new char[strlen(name)+1];
  r->hierarchy    = new char[strlen(hierarchy)+1];
  r->numerator    = new char[strlen(numerator)+1];
  r->denominator  = new char[strlen(denominator)+1];
  r->valid_phase_rate    = false;
  strcpy(r->name,name);
  strcpy(r->hierarchy,hierarchy);
  strcpy(r->numerator,numerator);
  strcpy(r->denominator,denominator);
  r->multiplier = multiplier;
  rate_map[name] = r;
}

void stats_t::register_phase_rate(const char* name, const char* hierarchy, const char* numerator, const char* denominator, double multiplier){
  // If the counter has been declared, mark it as a phase counter
  if(rate_map.find(name) != rate_map.end()){
    rate_map[name]->valid_phase_rate = true;
  }
  // If it does not exist, declare it and mark it as a phase counter
  else { 
    rate_t* r       = new rate_t;
    r->name         = new char[strlen(name)+1];
    r->hierarchy    = new char[strlen(hierarchy)+1];
    r->numerator    = new char[strlen(numerator)+1];
    r->denominator  = new char[strlen(denominator)+1];
    r->valid_phase_rate  = true;
    strcpy(r->name,name);
    strcpy(r->hierarchy,hierarchy);
    strcpy(r->numerator,numerator);
    strcpy(r->denominator,denominator);
    r->multiplier = multiplier;
    rate_map[name] = r;
  }
}


void stats_t::register_knob(const char* name, const char* hierarchy, unsigned int value){
  knob_t* k       = new knob_t;
  k->name         = new char[strlen(name)+1];
  k->hierarchy    = new char[strlen(hierarchy)+1];
  k->value        = value;
  strcpy(k->name,name);
  strcpy(k->hierarchy,hierarchy);
  knob_map[name] = k;
}


void stats_t::update_counter(const char* name,unsigned int inc){
  // If the counter has been declared and initialized
  if(counter_map.find(name) != counter_map.end()){
    counter_map[name]->count++;
    counter_map[name]->phase_count++;
  }
  // Tick the phase check mechanism if updating the 
  // counter on which phases are based on. Normally this
  // would be commit_count or cycle_count.
  if(!std::strcmp(name, phase_counter_name)){
    phase_tick();
  }
}

uint64_t stats_t::get_counter(const char* name){
  return counter_map[name]->count;
}

unsigned int stats_t::get_knob(const char* name){
  return knob_map[name]->value;
}

void stats_t::phase_tick(){
  if(counter_map[phase_counter_name]->phase_count >= phase_interval){
    phase_id++;
    update_rates();
    dump_phase_counters();
    dump_phase_rates();
    //dump_counters();
    //dump_rates();
    reset_phase_counters();
    fflush(0);
  }
}

void stats_t::update_rates(){
  std::map<std::string, rate_t*, ltstr>::iterator rate_iter;
  for(rate_iter = rate_map.begin();rate_iter != rate_map.end(); rate_iter++){
    if(counter_map[rate_iter->second->denominator]->count == 0){
      rate_iter->second->rate = (double)0.0;
    } else {
      rate_iter->second->rate = rate_iter->second->multiplier*
                                double(counter_map[rate_iter->second->numerator]->count)/
                                double(counter_map[rate_iter->second->denominator]->count);
    }

    if(counter_map[rate_iter->second->denominator]->phase_count == 0){
      rate_iter->second->phase_rate = (double)0.0;
    } else {
      rate_iter->second->phase_rate = rate_iter->second->multiplier*
                                      double(counter_map[rate_iter->second->numerator]->phase_count)/
                                      double(counter_map[rate_iter->second->denominator]->phase_count);
    }
  }
}

void stats_t::dump_counters(){
  fprintf(stats_log,"[stats]\n");
  std::map<std::string, counter_t*, ltstr>::iterator ctr_iter;
  for(ctr_iter = counter_map.begin();ctr_iter != counter_map.end(); ctr_iter++){
    fprintf(stats_log,"%s : %" PRIu64 "\n",ctr_iter->second->name, ctr_iter->second->count);
  }
}

void stats_t::dump_rates(){
  fprintf(stats_log,"[rates]\n");
  std::map<std::string, rate_t*, ltstr>::iterator rate_iter;
  for(rate_iter = rate_map.begin();rate_iter != rate_map.end(); rate_iter++){
    fprintf(stats_log,"%s : %2.2f\n",rate_iter->second->name, rate_iter->second->rate);
  }
}

void stats_t::dump_phase_counters(){
  fprintf(phase_log,"-------- Phase Counters Phase ID %" PRIu64 "--------\n",phase_id);
  std::map<std::string, counter_t*, ltstr>::iterator ctr_iter;
  for(ctr_iter = counter_map.begin();ctr_iter != counter_map.end(); ctr_iter++){
    if(ctr_iter->second->valid_phase_counter)
      fprintf(phase_log,"%s : %" PRIu64 "\n",ctr_iter->second->name, ctr_iter->second->phase_count);
  }
}

void stats_t::dump_phase_rates(){
  fprintf(phase_log,"-------- Phase Rates Phase ID %" PRIu64 "--------\n",phase_id);
  std::map<std::string, rate_t*, ltstr>::iterator rate_iter;
  for(rate_iter = rate_map.begin();rate_iter != rate_map.end(); rate_iter++){
    if(rate_iter->second->valid_phase_rate)
      fprintf(phase_log,"%s : %2.2f\n",rate_iter->second->name, rate_iter->second->phase_rate);
  }
}

void stats_t::dump_knobs(){
  fprintf(stats_log,"[knobs]\n");
  std::map<std::string, knob_t*, ltstr>::iterator knb_iter;
  for(knb_iter = knob_map.begin();knb_iter != knob_map.end(); knb_iter++){
    fprintf(stats_log,"%s : %u\n",knb_iter->second->name, knb_iter->second->value);
  }
}



