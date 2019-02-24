#ifndef SEQ
#define SEQ

//forward declaration
class system_unit_t;
class system_sim_t;

struct list_entry_t;

typedef struct{
  unsigned long timer;
  unsigned long max_strands;
  unsigned long num_channels;
  unsigned long base_timeout;
  unsigned long seq_dec_buffer_size;
  unsigned long prep_seq_buffer_size;
  unsigned long num_sequencers;
  struct list_entry_t* seq_dec_buffer;
  struct list_entry_t* prep_seq_buffer;
  system_sim_t* _system;
} sequencer_params_t; //bundled sequencer parameters


class sequencer_unit_t: public system_unit_t{
 public:
  unsigned long max_strands; //max number of reads for the sequencer 
  unsigned long used_strands; //number of strands used up on the sequencer 
  unsigned long wasted_strands; //junk strands that we are not interested in
  sequencer_unit_t(unsigned long max_sequencing, unsigned long num_channels);
};


//logical sequencer unit that envelopes all of the sequencers in the system
class sequencer_t{
 public:
  sequencer_t(sequencer_params_t sequencer_params);
  ~sequencer_t();
  void sequencer_stage(void);//top level interface to the system simulator
 private:
  sequencer_unit_t** sequencer_set; //set of sequencers in the system
  unsigned long num_sequencers; //number of sequencers in the system
  system_sim_t* _system; //pointer to the system simulator
  struct list_entry_t* seq_dec_buffer; //seq_dec_buffer pointer
  struct list_entry_t* prep_seq_buffer; // pointer to the prep_seq_buffer
  unsigned long seq_dec_buffer_size; //size of the seq_dec_buffer
  unsigned long prep_seq_buffer_size; //size of the prep_seq_buffer
  unsigned long base_timer; //time value used to refresh sequencer run times
  unsigned long base_timeout; //time value used to refresh the timeout timers

  void sequencer_backend(void);//this function will check sequencers and put transactions into the seq_dec_buffer
  void sequencer_frontend(void); //move transactions from the prep_seq_buffer to an open sequencer
  void sequencer_kickoff(void); //function that checks to see if sequencers are ready to be started
  void sequencer_timeoutstep(void); //timestep the the timeout counter
  void sequencer_submit(unsigned long transaction_ID, unsigned long sequencer_ID, unsigned long prep_seq_buffer_index); //submit a transaction to the sequencer
  int sequencer_avail(unsigned long transaction_ID); //check for an available sequencer
  void sequencer_complete(unsigned long sequencer_ID); //complete a finished sequencer
  int get_seqdec(unsigned long transaction_ID); //look for an open seq_dec_buffer entry
  void sequencer_timestep(unsigned long sequencer_ID); //decrement the sequencer's timer
};









#endif
