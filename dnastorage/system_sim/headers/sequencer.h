#ifndef SEQ
#define SEQ

//forward declaration
class system_unit;

class sequencer: public system_unit{
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
 sequencer(unsigned long timer,unsigned long max_sequencing, unsigned long num_channels, system* system_descriptor):system_unit(timer,num_channels,system_descriptor)

}



#endif
