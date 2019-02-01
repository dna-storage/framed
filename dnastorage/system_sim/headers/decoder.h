#ifndef DEC
#define DEC

class system_unit;

class decoder: public system_unit{
 public:
 decoder(unsigned long timer, unsigned long num_channels, system* system_descriptor) : system_unit(timer,num_channels,system_descriptor);
 
}

#endif
