#include"buffer.h"
#include<stdlib.h>

buffer_t::buffer_t(unsigned long buffer_size){
  this->buff=(list_entry_t*)malloc(sizeof(list_entry_t)*buffer_size);
  //initialize lists
  for(int i=0; i<buffer_size; i++){
    this->buff[i].used=0;
  }

  this->buffer_size=buffer_size;
}

buffer_t::~buffer_t(){
  free(this->buff);
}

void buffer_t::init_entry(unsigned long transaction_ID,unsigned long buffer_index){
  this->buff[buffer_index].transaction_index=transaction_ID;
  this->buff[buffer_index].used=1;
}

int buffer_t::get_free(void){
  
  for(int i=0; i<this->buffer_size; i++){
    if(this->buff[i].used==0) return i;
  }
  return -1; //failed to find a free spot

}

void buffer_t::operator++(){
  this->iterator++;
}

void buffer_t::iter_start(void){
  iterator=0;
}

list_entry_t* buffer_t::operator() (){
  if(this->iterator==this->buffer_size)return NULL;
  else return &(this->buff[this->iterator]);
}

void buffer_t::free_entry(unsigned long buffer_index){
  this->buff[buffer_index].used=0;
}

unsigned long buffer_t::get_iterator(void){
  return this->iterator;
}

int buffer_t::test_transaction(unsigned long transaction_ID, unsigned long& buffer_index){
  buffer_index=-1;
  for(int i=0; i<this->buffer_size;i++){
    if(this->buff[i].used==0) buffer_index=i;
    if(this->buff[i].used && this->buff[i].transaction_index==transaction_ID){
      buffer_index=i;
      return 1;
    }
  }
  return 0;
}
