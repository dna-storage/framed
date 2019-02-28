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
#include <stdlib.h>
#include <stdio.h>
#include <random>

generator_t::generator_t(generator_params_t generator_params)
{
  this->rate=generator_params.rate;
  this->max_file_size=generator_params.max_file_size;
  this->min_file_size=generator_params.min_file_size;
  this->unique_pools=generator_params.unique_pools;
  this->random_seed=generator_params.random_seed;
  this->_system=generator_params._system;
  gen=&generator_t::gen_poisson;
  //set up random number generators
  this->poisson_transactions = new std::poisson_distribution<unsigned long>(rate);
  this->rand_tool = new std::default_random_engine(this->random_seed);

  this->rand_pool= new std::uniform_int_distribution<unsigned long>(0,this->unique_pools-1);
  this->rand_file= new std::uniform_int_distribution<unsigned long>(this->min_file_size,this->max_file_size);
  
}

generator_t::~generator_t(){
  delete this->rand_tool;
  delete this->rand_file;
  delete this->rand_pool;
  delete this->poisson_transactions;
}


void generator_t::generator_stage(void){
  (this->*gen)();
}

//model transaction generation as a poisson process
void generator_t::gen_poisson(void){
  unsigned long number_transactions;
  trace_t* trace_transaction;
  std::default_random_engine def=*(this->rand_tool);
  std::poisson_distribution<unsigned long> _trans=*(this->poisson_transactions);
  std::uniform_int_distribution<unsigned long> _file=*(this->rand_file);
  std::uniform_int_distribution<unsigned long> _pool=*(this->rand_pool);
  
  //create some number of transactions based on the poisson process
  number_transactions=_trans(def);

  //create number_transactions trace_t objects
  for(int i=0; i<number_transactions; i++){
    //initialize the transactions
    trace_transaction=(trace_t*)malloc(sizeof(trace_t));
    trace_transaction->pool_ID=_pool(def);
    trace_transaction->time_stamp=_system->timer_tick;
    trace_transaction->file_size=_file(def);
    //add the transaction to the queue
    _system->queue_append(trace_transaction);
  }

}
