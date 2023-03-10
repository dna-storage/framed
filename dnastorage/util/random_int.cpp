#include <Python.h>
#include <time.h>
#include <random>

std::mt19937 RANDOM_NUMBER_SOURCE; 


//Python module seed C random number generator
static PyObject* generate_seed(PyObject *self){
  static std::random_device rd;
  return Py_BuildValue("I", rd()); //get a seed
}

static PyObject* set_seed(PyObject *self, PyObject* args){
  unsigned int seed;
  if(!PyArg_ParseTuple(args,"I",&seed)){
    return NULL;
  }
  RANDOM_NUMBER_SOURCE.seed(seed); //set the seed
  return Py_BuildValue("s",NULL);

}


//Python module to help with random integers, generates a random integer in a specified range
static PyObject* generate_rand_in_range(PyObject *self, PyObject *args){
  int lower,upper;
  if(!PyArg_ParseTuple(args,"ii",&lower,&upper)){
    return NULL;
  }
  std::uniform_int_distribution<int> dist(lower,upper);
  return Py_BuildValue("i",dist(RANDOM_NUMBER_SOURCE));
}

//Python module function to generate a random number between [0,1] (double)
static PyObject* generate_rand(PyObject *self){
  std::uniform_real_distribution<double> dist(0,1);
  return Py_BuildValue("d",dist(RANDOM_NUMBER_SOURCE));
}

static PyMethodDef generate_methods[] = {
   {"rand_in_range", (PyCFunction)generate_rand_in_range, METH_VARARGS, NULL },
   {"set_seed",(PyCFunction)set_seed,METH_VARARGS,NULL},
   {"seed",(PyCFunction)generate_seed,METH_NOARGS,NULL},
   {"rand",(PyCFunction)generate_rand,METH_NOARGS,NULL},
   { NULL, NULL, 0, NULL}
};

static struct PyModuleDef generate={
				    PyModuleDef_HEAD_INIT,
				    "generate",
				    NULL,
				    -1,
				    generate_methods
};

PyMODINIT_FUNC PyInit_generate(void){
  return PyModule_Create(&generate);
}
