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

  //FIX ME: Need this to change to allow different policies regarding the pools that can be in the sequencer at one time
  this->sequencer_poolpolicy=&sequencer_t::single_pool;
 
}



 
//look through the sequencers and find finished ones
//implements the backend completion of sequenced transactions
void sequencer_t::sequencer_backend(void){
  sequencer_unit_t* _sequencer;
  //iterate through the sequencers
  for(unsigned long i=0;i<(this->num_sequencers);i++){
    _sequencer=this->sequencer_set[i];
    unsigned long _timer=_sequencer->timer;
    int active=_sequencer->unit_active;
    if(active && _timer==0){
      //found a done sequencer, need to move each transaction in the list to the seq_dec_buffer
      this->sequencer_complete(i); //this function handles all of the steps needed for opening a decoder
    }
    else if(active && _timer>0){
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
  sequencer_unit_t* _sequencer;
  unsigned long transaction_strands; //number of strands that are still needed to be sequenced
  transaction_t* _window=this->_system->window;
  list_entry_t* _p_s=this->prep_seq_buffer;
  unsigned long transaction_ID;
  //iterate over all of the transactions in the prep_seq_buffer
  for(int i=0; i<QUEUE_SIZE; i++){
    if(_p_s[i].used){
      transaction_ID=_p_s.transaction_index;
      int sequencer_ID=this->sequencer_avail(transaction_ID);
      if(sequencer_ID>=0){
	this->sequencer_submit(transaction_ID,sequencer_ID); //submit the transaction to the sequencer
      }
    }
  }
  //kickoff sequencers
  this->sequencer_kickoff();
  //decrement standby timers
  this->sequencer_standbystep();
}

//find sequencers that are setup to be kicked off
void sequencer_t::sequencer_kickoff(void){
  sequencer_unit_t* _sequencer;

  //iterate through the sequencers and find expired stanby timers, or if the sequencing space has been exhausted
  for(int i=0; i<this->num_sequencers; i++){
    _sequencer=this->sequencer_set[i];
    if(_sequencer->used_strands==_sequencer->max_strands || (_sequencer->standby_timer==0 && _sequencer->next_open!=0)){
      _sequencer->unit_active=1; //activate the unit
      _sequencer->timer=base_timer-1; //initialize the timer for the sequencer
    }
  }
}

//decrease the standby timer 
void sequencer_t::standbystep(void){
  sequencer_unit_t* _sequencer;
  for(int i=0; i<this->num_sequencers; i++){
    _sequencer=this->sequencer_set[i];
    if(_sequencer->standby_timer>0 && !_sequencer->unit_active && _sequencer->next_open!=0){
      sequencer->standby_timer--;// decrement the stanby timer
    }
  }

}


//function to place a transaction into a sequencer
void sequencer_t::sequencer_submit(unsigned long transaction_ID, unsigned long sequencer_ID){
  	//found a sequencer that works with the transaction, need to place the transaction in that sequencer
	_sequencer=this->sequencer_set[sequencer_ID];
	transaction_strands=_window[transaction_ID].strands_to_sequence;
	if(transaction_strands>(_sequencer->max_strands-_sequencer->used_strands)){
	  //have more transaction strands than whats available in the sequencer
	  _sequencer->used_strands=_sequencer->max_strands;
	  _window[transaction_ID].strands_to_sequence-=(_sequencer->max_strands-_sequencer->used_strands);

	  //put the undesired strands in the sequencer 
	  if(_window[transaction_ID].undesired_strands_sequenced>0){
	    if(_window[transaction_ID].undesired_strands_sequenced>_sequencer->max_strands-_sequencer->used_strands){
	      _window[transaction_ID].undesired_strands-=(_sequencer->max_strands-_sequencer->used_strands);
	      _sequencer->wasted_strands+=(_sequencer->max_strands-_sequencer->used_strands);
	    }
	    else{
	      _sequencer->wasted_strands+=_window[transaction_ID].undesired_strands;
	      _window[transaction_ID].undesired_strands=0;
	    }
	  }
	  _window[transaction_ID].cracked++; //split transaction up and take note of it
	}
	else{
	  //transaction will fit in the remaining sequencer
	  _sequencer->used_strands+=_window[transaction_ID].strands_to_sequence;
	  _sequencer->wasted_strands+=_window[transaction_ID].undesired_strands_sequenced;
	  //free up the _p_s entry
	  _p_s[i].used=0;
	}
     	//make some changes to the standby timer and add the transactions to the sequencer
	_sequencer->transaction_slots[_sequencer->next_open]=transaction_ID;
	//if we are the first ones here, start up the standby timer
	if(_sequencer->next_open==0)_sequencer->standby_timer=this->base_standby_timer;
	_sequencer->next_open++;

}



//this function tries to find a sequencer that will work with transaction requested
//if there is no possible sequencer, then return -1
int sequencer_t::sequencer_avail(unsigned long transaction_ID){
  sequencer_unit_t* _sequencer;
  
  //iterate over all of the sequencers
  for(unsigned long i=0;i<(this->num_sequencers);i++){
    _sequencer=this->sequencer_set[i];
    //found an inactive sequencer
    if(!_sequencer->unit_active && _sequencer->used_strands < _sequencer->max_strands
       && _sequencer->next_open<_sequencer->num_channels){
      if(_sequencer->next_open==0){
	//know that this sequencer is unused and unscheduled by any other transaction
	return i;
      }
      else{
	//need to look at other transactions and based off of policy all the requested transaction in
 	if(this->sequencer_poolpolicy(transaction_ID,i)) return i;
      }
    }

  }
}


//policy function for the sequencer, only lets 1 pool be in a single sequencer
int sequencer_t::single_pool(unsigned long transaction_ID, unsigned long sequencer_ID){
  transaction_t* _window=this->_system->window;
  sequencer_unit_t _sequencer=this->sequencer_set[sequencer_ID];
  unsigned long sequencer_trans;
  
  //iterate through all the transactions on a given sequencer and evaluate the pools
  for(int i=0; i< _sequencer->next_open ; i++){
    sequencer_trans=_sequencer->transaction_slots[i];
    if(_window[sequencer_trans].pool_ID==_window[transaction_ID].pool_ID) return 0;
  }
  return 1;

}



//complete the decoder by placing the transactions into the seq_dec_buffer
void sequencer_t::sequencer_complete(unsigned long sequencer_ID){
  unsigned long number_transactions=(this->sequencer_set[sequencer_ID])->next_open;
  sequencer_unit_t* _sequencer=this->sequencer_set[sequencer_ID];
  list_entry_t* _s_d=this->seq_dec_buffer;
  transaction_t* _window=this->_system->window;
  
  for(_sequencer->transaction_pointer;
      _sequencer->transaction_pointer < number_transactions;
      _sequencer->transaction_pointer++){
    
    unsigned long i=_sequencer->transaction_pointer;
    unsigned long transaction_number=_sequencer->transaction_slots[i];
    int seq_dec_index= this->get_seqdec(transaction_number);
    if(seq_dec_index==-1) return; //no more buffer spots 
    
    if(seq_dec_index>=0){ //need a new buffer spot
      //put the transactions that can fit into the seq_dec_buffer, and mark it used
      _s_d[seq_dec_index].used=1;
      _s_d[seq_dex_index].transaction_index=transaction_number;
      _window[transaction_number].cracked_count--;
    }
    else{
      //transaction already in the buffer, just decrement the cracked count for the transaction
      _window[transaction_number].cracked_count--;
    }
    
  }
  //if we made it here, we moved all of the transactions from the sequencer to the buffer
  _sequencer->unit_active=0; //kill the activity of the unit
  _sequencer->used_strands=0; //reset the used strands
  _sequencer->next_open=0; //reset the next open pointer
  _sequencer->transaction_pointer=0;//reset the transaction pointer
}




//look through the seqdec buffer for an open spot.
// --open spot required if a transaction already exists in the seq_dec_buffer. This suggests a previously cracked transaction
// --Return -1 when an open spot in the buffer cannot be found, this will halt the sequencer from becoming unactivated
int sequencer_t::get_seqdec(unsigned long transaction_ID){
  list_entry_t* _s_d=this->seq_dec_buffer;
  
  
  //make syre that the seq_dec_buffer does not already have the transaction, do not want to double place transactions
  for(int i=0; i<QUEUE_SIZE;i++){
    if(transaction_ID==_s_d[i].transaction_index && _s_d[i].used) return -9999;
  }
  //find an open seq_dec_buffer location
  for(int i=0; i<QUEUE_SIZE;i++){
    if(_s_d[i].used==0) return i; // return an unused list_entry
  }
  return -1;

}


//decreases the timer for a sequencer that is active
void sequencer_t::sequencer_timestep(unsigned long sequencer_ID){
  sequencer_unit_t* _sequencer=this->sequencer_set[sequencer_ID];
  _sequencer->timer--;
}
