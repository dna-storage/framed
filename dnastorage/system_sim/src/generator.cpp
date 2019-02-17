#include <stdlib.h>
#include <stdio.h>
#include <random>
#include "dna_storage_attributes.h"
#include "decoder.h"
#include "prep.h"
#include "sequencer.h"
#include "stdio.h"
#include "utlist.h"

generator_t::generator_t(float rate,unsigned long max_file_size,unsigned long min_file_size,
			 unsigned long unique_pools, int random_seed, system_sim_t* _system)
{
  this->rate=rate;
  this->max_file_size=max_file_size;
  this->min_file_size=min_file_size;
  this->unique_pools=unique_pools;
  this->random_seed=random_seed;
  this->_system=_system;
  srand(this->random_seed);//set the random generator seed
  gen=&generator_t::gen_poisson;
  
  //set up random number generators
  this->poisson_transactions = new std::poisson_distribution<int>(rate);
  this->rand_pois = new std::default_random_engine();
  this->rand_file = new std::default_random_engine();
  this->rand_pool = new std::default_random_engine();
  this->rand_file.min(min_file_size);
  this->rand_file.max(max_file_size);
  this->rand_pool.min(0);
  this->rand_pool.max(this->unique_pools);
}

generator_t::~generator_t(){
  delete this->rand_pois;
  delete this->rand_file;
  delete this->poisson_transactions;
}


void generator_t::generator_stage(void){
  this->gen();
}

//model transaction generation as a poisson process
void generator_t::gen_poisson(void){
  int number_transactions;
  trace_t* trace_transaction;
  //create some number of transactions based on the poisson process
  number_transactions=this->poisson_transactions(this->rand_pois);

  //create number_transactions trace_t objects
  for(int i=0; i<number_transactions; i++){
    //initialize the transactions
    trace_transaction=(trace_t*)malloc(sizeof(trace_t));
    trace_transaction->pool_ID=this->rand_pool();
    trace_transaction->time_stamp=_system->timer_tick;
    trace_transaction->file_size=this->rand_file();
    //add the transaction to the queue
    _system->queue_append(trace_transaction);
  }

}
