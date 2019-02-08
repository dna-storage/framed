#include "dna_storage_attributes.h"
#include "decoder.h"
#include "prep.h"
#include "sequencer.h"
#include <stdlib.h>
sequencer_unit_t::sequencer_unit_t(unsigned long timer,unsigned long max_sequencing,
				   unsigned long num_channels):system_unit_t(num_channels){
  this->max_strands=max_sequencing;
  this->used_strands=0;
  this->wasted_strands=0;
  this->utilization=0;
}


sequencer_t:: sequencer_t(unsigned long timer,unsigned long max_strands,
			  unsigned long num_channels,unsigned long buffer_size,
			  unsigned long num_sequencers,list_entry_t* seq_dec_buffer,
			  list_entry_t* prep_seq_buffer, system_sim_t* _system,
			  unsigned long base_standby_timer){

  
  this->num_sequencers=num_sequencers;
  this->_system=_system;
  this->seq_dec_buffer=seq_dec_buffer;
  this->prep_seq_buffer=prep_seq_buffer;
  this->buffer_size=buffer_size;
  this->base_timer=timer;
  this->base_standby_timer=base_standby_timer;
  //allocate sequencer units
  this->sequencer_set=(sequencer_unit_t**)malloc(sizeof(sequencer_unit_t*)*this->num_sequencers);
  for(int i=0; i<this->num_sequencers; i++) this->sequencer_set[i]=new sequencer_unit_t(max_strands,0);
 
}



 
//look through the sequencers and find finished ones
//implements the backend completion of sequenced transactions
void sequencer_t::sequencer_backend(void){
  for(unsigned long i=0;i<(this->num_sequencers);i++){
    unsigned long timer=(this->sequencer_set[i])->timer;
    int active=(this->sequencer_set[i])->unit_active;
    if(active && timer==0){
      //found a done sequencer, need to move each transaction in the list to the seq_dec_buffer
      this->sequencer_complete(i); //this function handles all of the steps needed for opening a decoder
    }
    else if(active && timer>0){
      //active sequencer, decrement the timer
      this->sequencer_timestep(i);
    }

  }
}

//this function implements the frontend of the sequencer
//the frontend is responsible for moving transactions from the prep_seq_buffer to a sequencing station
/*
  Steps:
  1. Schedule to a sequencer. If the sequencer was totally idle (no other transaction scheduled yet) start the standby timer.
     If a sequencer has no been kicked off yet, but some other transaction has been scheduled to it then make sure the sequencer placement policy
     allows for it. If a transaction cannot fit in its entirety, need to track that a crack was made.
  2. If the standby counter is 0 for a sequencer, tell it to become active
  3. decrement standby counters
*/
void sequencer_t::sequencer_frontend(void){

  //iterate over all of the transactions in the prep_seq_buffer
  for(int i=0; i<QUEUE_SIZE; i++){
    if(this->prep_seq_buffer[i].used){
      unsigned long sequencer_ID=this->sequencer_avail(this->prep_seq_buffer[i].transaction_index);
      if(sequencer_ID>=0){
	//found a sequencer that works with the transaction 
	
      }
    }
  }
}

//this function tries to find a sequencer that will work with transaction requested
//if there is no possible sequencer, then return -1
int sequencer_t::sequencer_avail(unsigned long transaction_ID){
  //iterate over all of the sequencers
  for(unsigned long i=0;i<(this->num_sequencers);i++){
    //found an inactive sequencer
    if(!this->sequencer_set[i]->unit_active && this->sequencer_set[i]->used_strands<this->sequencer_set[i]->max_strands){
      if(this->sequencer_set[i]->next_open==0){
	//know that this sequencer is unused and unscheduled by any other transaction
	return i;
      }
      else{
	//need to look at other transactions and based off of policy all the requested transaction in
	if(this->sequencer_policy(transaction_ID,i)) return i;
      }
    }

  }


}



//complete the decoder by placing the transactions into the seq_dec_buffer
void sequencer_t::sequencer_complete(unsigned long sequencer_ID){
  unsigned long number_transactions=(this->sequencer_set[sequencer_ID])->next_open;
  
  for((this->sequencer_set[sequencer_ID])->transaction_pointer;
      (this->sequencer_set[sequencer_ID])->transaction_pointer<number_transactions;
      (this->sequencer_set[sequencer_ID])->transaction_pointer++){
    
    unsigned long i=(this->sequencer_set[sequencer_ID])->transaction_pointer;
    unsigned long transaction_number=transaction_slots[i];
    int seq_dec_index=(this->get_seqdec(transaction_number));
    if(seq_dec_index==-1) return; //no more buffer spots 
    
    if(seq_dec_index>=0){ //need a new buffer spot
      //put the transactions that can fit into the seq_dec_buffer, and mark it used
      this->seq_dec_buffer[seq_dec_index].used=1;
      this->seq_dex_buffer[seq_dex_index].transaction_index=transaction_number;
      (this->_system)->window[transaction_number].cracked_count--;
    }
    else{
      //transaction already in the buffer, just decrement the cracked count for the transaction
      (this->_system)->window[transaction_number].cracked_count--;
    }
    
  }
  //if we made it here, we moved all of the transactions from the sequencer to the buffer
  (this->sequencer_set[sequencer_ID])->unit_active=0; //kill the activity of the unit
  (this->sequencer_set[sequencer_ID])->used_strands=0; //reset the used strands
  (this->sequencer_set[sequencers_ID])->next_open=0; //reset the next open pointer
}




//look through the seqdec buffer for an open spot.
// --open spot required if a transaction already exists in the seq_dec_buffer. This suggests a previously cracked transaction
// --Return -1 when an open spot in the buffer cannot be found, this will halt the sequencer from becoming unactivated
int sequencer_t::get_seqdec(unsigned long transaction_ID){

  //make syre that the seq_dec_buffer does not already have the transaction, do not want to double place transactions
  for(int i=0; i<QUEUE_SIZE;i++){
    if(transaction_ID==(this->seq_dec_buffer.transaction_index) && this->seq_dec_buffer[i].used) return -9999;
  }
  //find an open seq_dec_buffer location
  for(int i=0; i<QUEUE_SIZE;i++){
    if((this->seq_dec_buffer[i].used==0)) return i; // return an unused list_entry
  }
  return -1;

}



void sequencer_t::sequencer_timestep(unsigned long sequencer_ID){
  (this->sequencer_set[sequencer_ID])->timer--;
}
