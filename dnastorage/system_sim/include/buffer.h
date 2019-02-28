#ifndef BUFF_H
#define BUFF_H

typedef struct list_entry_t{
  unsigned long transaction_index; //index into the window[] 
  int used; //indicates whether the entry is used or not
} list_entry_t; //structure that defines each buffer's entry 


class buffer_t{ //class that implements support for connecting buffers
  public:
  list_entry_t* buff;
  unsigned long buffer_size;

  
  buffer_t(unsigned long);
  ~buffer_t();
  void operator++(void);
  list_entry_t* operator()(void);
  void iter_start(void);
  void init_entry(unsigned long transaction_ID, unsigned long buff_index);
  int get_free(void);
  unsigned long get_iterator(void);
  void free_entry(unsigned long buffer_index);
  int test_transaction(unsigned long transaction_ID, unsigned long& buffer_index);
 private:
  unsigned long iterator;
      
};



#endif 
