#ifndef PREP
#define PREP

class system_unit;

class prep: public system_unit{
 public:
 prep(unsigned long timer, unsigned long num_channels, system* system_descriptor) : system_unit(timer,num_channels,system_descriptor)
 
}

#endif
