#ifndef STORAGE
#define STORAGE



typedef struct{
  unsigned long timer; //timer keeping track of whether the pool is usable
  unsigned long remaining_reads; //keeps track of the remaining reads of this pool
  unsigned long in_use; //flag that indicates if pool can be used, make sure back to back uses are blocked
} pool_char_t; //structure that defines the characteristics of each copy

//structure used to model a pool containing files
typedef struct{
  pool_char_t* copies; //array of copies
} pool_model_t; //structure to model a pool in the storage system




typedef struct{
  float sequencing_efficiency;
  unsigned long average_pool_capacity;
  unsigned long number_pools;
  float  bytes_per_strand;
  unsigned long pool_copies;
  unsigned long pool_write_time;
  unsigned long pool_wait_time;
  unsigned long number_reads;
} storage_params_t; //bundled storage parameters




//class that models the system storage
class system_storage_t{
 public:

   system_storage_t(storage_params_t storage_params);
  ~system_storage_t();
  int storage_poolavailable(unsigned long pool_ID);
  void storage_readmanage(unsigned long pool_ID, unsigned long copy_ID);
  void storage_poolrestore(void);

 private:    
  unsigned long number_pools; //number of unique pools in the storage system 
  unsigned long average_pool_capacity; //average capacity of the pool (MB), need to model whole pool sequencing 
  unsigned long pool_copies;//number of copies each pool has
  unsigned long number_reads;//number of reads that a pool can withstand 
  unsigned long pool_write_time;//time it takes to replenish pool after using all of its reads
  unsigned long pool_wait_time;//time it takes for a pool to be available after using it
  pool_model_t* pools; //array of unique pools in the system
  
};







#endif 
