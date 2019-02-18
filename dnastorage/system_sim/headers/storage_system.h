#ifdefine STORAGE
#define STORAGE
struct pool_model_t;

typedef struct{
  float sequencing_efficiency;
  unsigned long average_pool_capacity;
  unsigned long number_pools;
  unsigned long bytes_per_strand;
  unsigned long pool_copies;
  unsigned long pool_write_time;
  unsigned long pool_wait_time;
} storage_params_t;




//class that models the system storage
class system_storage_t{
 public:
  unsigned long number_pools;
  unsigned long average_pool_capacity; //average capacity of the pool (MB), need to model whole pool sequencing 
  float sequencing_efficiency; //models how well a certain 
  //parameters of pool re-use
  unsigned long pool_copies;//number of copies the pool has
  unsigned long number_reads;//number of reads that a pool can withstand 
  unsigned long pool_write_time;//time it takes to replenish pool after using all of its reads
  unsigned long pool_wait_time;//time it takes for a pool to be available after using it
  //pool object that holds files
  struct pool_model_t* pools;
  system_storage_t(storage_params_t storage_params);
  ~system_storage_t();
  int storage_poolavailable(unsigned long pool_ID);
  void storage_readmanage(unsigned long pool_ID, unsigned long copy_ID);
  void storage_poolrestore(void);
};







#endif 
