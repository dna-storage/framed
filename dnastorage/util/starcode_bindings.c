#include <Python.h>
#include <assert.h>
#include <stdlib.h>
#include "starcode.h"


//Python module to help with random integers, generates a random integer in a specified range
static PyObject* starcode_top(PyObject *self, PyObject *args){
  int ret_status=0;
  int number_strands; //number of strands that are part of the clustering set
  int clusteralg; //Sphere clustering/message passing clustering/...
  PyObject* strand_array; //pointer to python object that contains the input aray data set
  char** input_strands; //array of strands that need to be clustered
  int tau; //Max levenshtein distance
  double parent_to_child; //ratio of parent to child reads to consider the parent as the canonical strand
  int thrmax; //number of threads to run starcode with

  
  if(!PyArg_ParseTuple(args,"iiidOi",&tau,&thrmax,&clusteralg,&parent_to_child,&strand_array,&number_strands)){
    return NULL;
  }
  
  input_strands=(char**)malloc(sizeof(char*)*number_strands);
  assert(input_strands!=NULL);
  assert(tau<STARCODE_MAX_TAU);//make sure max LD is in range

  //need to take the strand list python object and put it into a correct C type form
  for(int i=0;i<number_strands;i++){
    PyObject* item=PySequence_Fast_GET_ITEM(strand_array,i);
    if(!item){
      Py_XDECREF(strand_array);
      free(input_strands);
      assert(0);
    }
    input_strands[i]=(char*)malloc(strlen(PyString_AsString(item))+1); //allocate space for a c string
    strcpy(input_strands[i],PyString_AsString(item)); //copy string from Python's memory space
    //Py_XDECREF(item);
  }
  //Py_XDECREF(strand_array);
 
  
  //call the starcode algorithm
  ret_status=starcode(NULL,NULL,NULL,NULL,tau,1,thrmax,clusteralg,parent_to_child,0,0,0,input_strands,number_strands);

  assert(ret_status==0); //make sure the algorithm finished successfully
  for(int i=0; i<number_strands;i++) free(input_strands[i]);
  free(input_strands);
  return Py_BuildValue("i",ret_status);
  
 
}


//Python module function to generate a random number between [0,1] (double)
static PyObject* starcode_get_cluster(PyObject *self){
  //gets a cluster from starcode
  char** return_cluster;

  //printf("Get cluster from starcode\n");
  return_cluster=get_cluster();
  //printf("retrieved cluster from starcode\n");
  if(return_cluster==NULL){
    //need to signal that there are no more clusters left to look at
    Py_INCREF(Py_None);
    return Py_None;
  }
  else{
    //need to create a python list object from the array of strings
    int i=0;
    PyObject* string_list=PyList_New(0);
    while(return_cluster[i]!=NULL){
      //printf("strand in cluster: %s\n",return_cluster[i]);
      PyObject* string_item=Py_BuildValue("s",return_cluster[i]);
      int append_results=PyList_Append(string_list,string_item);
      assert(!append_results);
      i++;
    }
    return string_list;
  }
 
}


static PyMethodDef starcode_methods[] = {
   { "starcode_top", (PyCFunction)starcode_top, METH_VARARGS, NULL },
   {"starcode_get_cluster",(PyCFunction)starcode_get_cluster,METH_NOARGS,NULL},
   { NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC initstarcode_bindings() {
   Py_InitModule3("starcode_bindings", starcode_methods, "Starcode clustering algorithm interface");
}
