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
  unsigned long buffer_size;
  unsigned long num_sequencers;
  struct list_entry_t* seq_dec_buffer;
  struct list_entry_t* prep_seq_buffer;
  system_sim_t* _system;
} sequencer_params_t;


class sequencer_unit_t: public system_unit_t{
 public:
  //max number of strands the sequencer can handle
  unsigned long max_strands;
  //number of strands currently used on the current sequencing runs 
  unsigned long used_strands;
  //strands that are wasted, function of policy of sequencing (e.g sequence whole pool/sequencing efficiency)
  unsigned long wasted_strands;

  
  //measure of each sequencers utilization at a given time step (used_strands-wasted)/max_strands
  float utilization;

  //constructor for base class and sequencer class
  sequencer_unit_t(unsigned long max_sequencing, unsigned long num_channels);
};


//logical sequencer unit that envelopes all of the sequencers in the system
class sequencer_t{

 public:
  sequencer_unit_t** sequencer_set;
  unsigned long num_sequencers;
  system_sim_t* _system;
  struct list_entry_t* seq_dec_buffer;
  struct list_entry_t* prep_seq_buffer;
  unsigned long buffer_size;
  unsigned long base_timer;
  unsigned long base_timeout;
  sequencer_t(sequencer_params_t sequencer_params);
  ~sequencer_t();

  

  //top level sequencer functions to be called by the top level simulator
  void sequencer_backend(void);//this function will check sequencers and put transactions into the seq_dec_buffer
  void sequencer_frontend(void);

};









#endif
