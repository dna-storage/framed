#include "sequencer.h"

sequencer::sequencer(unsigned long timer,unsigned long max_sequencing, unsigned long num_channels, system* system_descriptor){
  this->max_strands=max_sequencing;
  this->used_strands=0;
  this->wasted_strands=0;
  this->utilization=0;
}
