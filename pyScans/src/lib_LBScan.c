#include "/usr/include/python2.7/Python.h"
//#include "/usr/share/pyshared/numpy/core/include/numpy/arrayobject.h"
//#include "/usr/share/pyshared/numpy/core/include/numpy/arrayscalars.h"
#include "/usr/lib64/python2.7/site-packages/numpy/core/include/numpy/arrayobject.h"
#include "/usr/lib64/python2.7/site-packages/numpy/core/include/numpy/arrayscalars.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//##############################################################################################
double **ptrvector(int n)  {
  double **v;
  v=(double **)malloc((size_t) (n*sizeof(double)));
  if (!v)   {
    printf("In **ptrvector. Allocation of memory for double array failed.");
    exit(0);  }
  return v;
}

PyArrayObject *pymatrix(PyObject *objin)  {
  return (PyArrayObject *) PyArray_ContiguousFromObject(objin,
							NPY_DOUBLE, 2,2);
}

double *pyvector_to_Carrayptrs(PyArrayObject *arrayin){
  int n;
  n=arrayin->dimensions[0];
  return (double *) arrayin->data;  /* pointer to arrayin data as double */
}

int *pyvector_to_Carrayptrs_int(PyArrayObject *arrayin){
  int n;
  n=arrayin->dimensions[0];
  return (int *) arrayin->data;  /* pointer to arrayin data as double */
}

double **pymatrix_to_Carrayptrs(PyArrayObject *arrayin)  {
  double **c, *a;
  int i,n,m;
  
  n=arrayin->dimensions[0];
  m=arrayin->dimensions[1];
  c=ptrvector(n);
  a=(double *) arrayin->data;  /* pointer to arrayin data as double */
  for ( i=0; i<n; i++)  {
    c[i]=a+i*m;  }
  return c;
}

void free_Carrayptrs(double *v){
  free((char*) v);
}

void free_Cint2Darrayptrs(double **v)  {
  free((char*) v);
}

int not_doublevector(PyArrayObject *vec){
  if (vec->descr->type_num != NPY_DOUBLE || vec->nd != 1)  {
    PyErr_SetString(PyExc_ValueError,
		    "In not_doublevector: array must be of type Float and 1 dimensional (n).");
    return 1;  }
  return 0;
}

//##############################################################################################
void matrix2x2_multi_xy( double *x,  double *y,  double phi)
{
  double tmp;
  tmp =  cos(phi) * (*x) + sin(phi) * (*y);
  *y = -sin(phi) * (*x) + cos(phi) * (*y);
  *x = tmp;
}

void matrix2x2_multi_xz( double *x,  double *z,  double theta)
{
  double tmp;
  tmp =  cos(theta)* (*x) + sin(theta)* (*z);
  *z = -sin(theta)* (*x) + cos(theta)* (*z);
  *x = tmp;
}
//##############################################################################################
static PyObject *LB_rotmatrix_multi(PyObject *self, PyObject *args)
{
  /* basic variables */
  long i;
  int num, dims[2];
  double pi=M_PI;
  //  double pi=4.*atan(1.), radeg=(180./pi);
  
  /* variables for inputs */
  PyArrayObject *theta_asp, *phi_asp; // theta and phi of anti-sun vector
  PyArrayObject *phi_antisun, *phi_boresight; // precession and spin axes of satellite rotation
  double theta_antisun, theta_boresight;

  /* variable for outputs */
  PyArrayObject *pout;
  
  /* Parse tuples separately since args will differ between C fcns */
  if (!PyArg_ParseTuple(args, "O!O!dO!dO!", 
			&PyArray_Type, &theta_asp, &PyArray_Type, &phi_asp,
			&theta_antisun, &PyArray_Type, &phi_antisun,
			&theta_boresight, &PyArray_Type, &phi_boresight))  return NULL;
  if (NULL == theta_asp)  return NULL;
  if (NULL == phi_asp)  return NULL;
  if (NULL == phi_antisun)  return NULL;
  if (NULL == phi_boresight)  return NULL;
  
  /* Check that objects are 'double' type and vectors
     Not needed if python wrapper function checks before call to this routine */
  if (not_doublevector(theta_asp)) return NULL;
  if (not_doublevector(phi_asp)) return NULL;
  if (not_doublevector(phi_antisun)) return NULL;
  if (not_doublevector(phi_boresight)) return NULL;
  
  /* Get vector dimension. */
  dims[0] = 3;
  num=dims[1]=theta_asp->dimensions[0];
  pout=(PyArrayObject *) PyArray_FromDims(2,dims,NPY_DOUBLE);

  /* Change contiguous arrays into C * arrays   */
  double *c_theta_asp, *c_phi_asp;
  double *c_phi_antisun, *c_phi_boresight; // precession axis, spin axis
  double **c_pout;
  c_theta_asp = pyvector_to_Carrayptrs(theta_asp);
  c_phi_asp = pyvector_to_Carrayptrs(phi_asp);

  c_phi_antisun = pyvector_to_Carrayptrs(phi_antisun);
  c_phi_boresight = pyvector_to_Carrayptrs(phi_boresight);
  c_pout = pymatrix_to_Carrayptrs(pout);

  //  printf("     > %lf %lf \n  ", theta_antisun/pi*180., theta_boresight/pi*180.);
  double xi, yi, zi, rel_phi;
  for(i=0;i<num;i++){
    
    xi = sin(c_theta_asp[i])*cos(c_phi_asp[i]);
    yi = sin(c_theta_asp[i])*sin(c_phi_asp[i]);
    zi = cos(c_theta_asp[i]);
    
    if (yi >= 0) rel_phi = acos(xi);
    if (yi < 0) rel_phi = -acos(xi) + 2.*pi;
    matrix2x2_multi_xy( &xi, &yi,  rel_phi);
    matrix2x2_multi_xz( &xi, &zi,  pi/2.);
    matrix2x2_multi_xz( &xi, &zi, -theta_boresight);
    matrix2x2_multi_xy( &xi, &yi, -c_phi_boresight[i]);
    matrix2x2_multi_xz( &xi, &zi, -theta_antisun);
    matrix2x2_multi_xy( &xi, &yi, -c_phi_antisun[i]);
    matrix2x2_multi_xz( &xi, &zi, -pi/2.);
    matrix2x2_multi_xy( &xi, &yi, -rel_phi);
    
    c_pout[0][i] = xi;
    c_pout[1][i] = yi;
    c_pout[2][i] = zi;
  }
  
//  free_Carrayptrs(c_theta_asp);
//  free_Carrayptrs(c_phi_asp);
//  free_Carrayptrs(c_phi_antisun);
//  free_Carrayptrs(c_phi_boresight);
  //  free_Carrayptrs(c_pout);
  //  printf("=> %lf %lf %lf \n", c_pout[0][200], c_pout[1][200], c_pout[2][200]);
  //  return Py_BuildValue("i", 1);
  return PyArray_Return(pout);
}
//##############################################################################################
static PyObject *LB_rotmatrix_multi_ellip(PyObject *self, PyObject *args)
{
  /* basic variables */
  long i;
  int num, dims[2];
  double pi=M_PI;
  //  double pi=4.*atan(1.), radeg=(180./pi);
  
  /* variables for inputs */
  PyArrayObject *theta_asp, *phi_asp; // theta and phi of anti-sun vector
  PyArrayObject *phi_antisun, *phi_boresight; // precession and spin axes of satellite rotation
  PyArrayObject *theta_antisun; 
  double theta_boresight;

  /* variable for outputs */
  PyArrayObject *pout;
  
  /* Parse tuples separately since args will differ between C fcns */
  if (!PyArg_ParseTuple(args, "O!O!O!O!dO!", 
			&PyArray_Type, &theta_asp, &PyArray_Type, &phi_asp,
			&PyArray_Type, &theta_antisun, &PyArray_Type, &phi_antisun,
			&theta_boresight, &PyArray_Type, &phi_boresight))  return NULL;
  if (NULL == theta_asp)  return NULL;
  if (NULL == phi_asp)  return NULL;
  if (NULL == phi_antisun)  return NULL;
  if (NULL == phi_boresight)  return NULL;
  
  /* Check that objects are 'double' type and vectors
     Not needed if python wrapper function checks before call to this routine */
  if (not_doublevector(theta_asp)) return NULL;
  if (not_doublevector(phi_asp)) return NULL;
  if (not_doublevector(phi_antisun)) return NULL;
  if (not_doublevector(phi_boresight)) return NULL;
  
  /* Get vector dimension. */
  dims[0] = 3;
  num=dims[1]=theta_asp->dimensions[0];
  pout=(PyArrayObject *) PyArray_FromDims(2,dims,NPY_DOUBLE);

  /* Change contiguous arrays into C * arrays   */
  double *c_theta_asp, *c_phi_asp;
  double *c_phi_antisun, *c_phi_boresight; // precession axis, spin axis
  double *c_theta_antisun;
  double **c_pout;
  c_theta_asp = pyvector_to_Carrayptrs(theta_asp);
  c_phi_asp = pyvector_to_Carrayptrs(phi_asp);

  c_phi_antisun = pyvector_to_Carrayptrs(phi_antisun);
  c_phi_boresight = pyvector_to_Carrayptrs(phi_boresight);
  c_theta_antisun = pyvector_to_Carrayptrs(theta_antisun);
  c_pout = pymatrix_to_Carrayptrs(pout);

  //  printf("     > %lf %lf \n  ", theta_antisun/pi*180., theta_boresight/pi*180.);
  double xi, yi, zi, rel_phi;
  for(i=0;i<num;i++){
    
    xi = sin(c_theta_asp[i])*cos(c_phi_asp[i]);
    yi = sin(c_theta_asp[i])*sin(c_phi_asp[i]);
    zi = cos(c_theta_asp[i]);
    
    if (yi >= 0) rel_phi = acos(xi);
    if (yi < 0) rel_phi = -acos(xi) + 2.*pi;
    matrix2x2_multi_xy( &xi, &yi,  rel_phi);
    matrix2x2_multi_xz( &xi, &zi,  pi/2.);
    matrix2x2_multi_xz( &xi, &zi, -theta_boresight);
    matrix2x2_multi_xy( &xi, &yi, -c_phi_boresight[i]);
    matrix2x2_multi_xz( &xi, &zi, -c_theta_antisun[i]);
    matrix2x2_multi_xy( &xi, &yi, -c_phi_antisun[i]);
    matrix2x2_multi_xz( &xi, &zi, -pi/2.);
    matrix2x2_multi_xy( &xi, &yi, -rel_phi);
    
    c_pout[0][i] = xi;
    c_pout[1][i] = yi;
    c_pout[2][i] = zi;
  }
  
//  free_Carrayptrs(c_theta_asp);
//  free_Carrayptrs(c_phi_asp);
//  free_Carrayptrs(c_phi_antisun);
//  free_Carrayptrs(c_phi_boresight);
  //  free_Carrayptrs(c_pout);
  //  printf("=> %lf %lf %lf \n", c_pout[0][200], c_pout[1][200], c_pout[2][200]);
  //  return Py_BuildValue("i", 1);
  return PyArray_Return(pout);
}
//##############################################################################################

static PyObject *Maps_summingup(PyObject *self, PyObject *args)
{
  /* basic variables */
  int i;
  int dims[2];
  
  /* variables for inputs */
  int nbPix;
  PyArrayObject *ipix, *alpha; // theta and phi of anti-sun vector
  
  /* variable for outputs */
  PyArrayObject *mapout;
  
  /* Parse tuples separately since args will differ between C fcns */
  if (!PyArg_ParseTuple(args, "iO!O!", 
			&nbPix, 
			&PyArray_Type, &ipix, &PyArray_Type, &alpha))  return NULL;
  if (NULL == ipix)  return NULL;
  if (NULL == alpha)  return NULL;
  if (not_doublevector(ipix)) return NULL;
  if (not_doublevector(alpha)) return NULL;
  
  /* Get vector dimension. */
  dims[0] = 7;
  dims[1] = nbPix;
  mapout=(PyArrayObject *) PyArray_FromDims(2,dims,NPY_DOUBLE);
  int nbTOD = alpha->dimensions[0]; 

  /* Change contiguous arrays into C * arrays   */
  double *c_ipix, *c_alpha;
  double **c_mapout;

  c_ipix = pyvector_to_Carrayptrs(ipix);
  c_alpha = pyvector_to_Carrayptrs(alpha);
  c_mapout = pymatrix_to_Carrayptrs(mapout);

  for(i=0;i<nbTOD;i++){
    c_mapout[0][(int)c_ipix[i]] += 1.;
    c_mapout[1][(int)c_ipix[i]] += cos(1.*c_alpha[i]);
    c_mapout[2][(int)c_ipix[i]] += sin(1.*c_alpha[i]);
    c_mapout[3][(int)c_ipix[i]] += cos(2.*c_alpha[i]);
    c_mapout[4][(int)c_ipix[i]] += sin(2.*c_alpha[i]);
    c_mapout[5][(int)c_ipix[i]] += cos(4.*c_alpha[i]);
    c_mapout[6][(int)c_ipix[i]] += sin(4.*c_alpha[i]);
  }
  
  //  printf("ok5\n");
  //    free((char*) c_alpha);
  //    free((char*) c_ipix);
  //  free((char**) c_mapout);
  //  free_Carrayptrs(c_alpha);
  //  free_Carrayptrs(c_ipix);
  //  free_Cint2Darrayptrs(c_mapout);
  //  printf("=> %lf %lf %lf \n", c_pout[0][200], c_pout[1][200], c_pout[2][200]);
  //  return Py_BuildValue("i", 1);
  return PyArray_Return(mapout);
}

//##############################################################################################

static PyMethodDef methods[] = {
  {"LB_rotmatrix_multi", LB_rotmatrix_multi, METH_VARARGS, ""},
  {"LB_rotmatrix_multi_ellip", LB_rotmatrix_multi_ellip, METH_VARARGS, ""},
  {"Maps_summingup", Maps_summingup, METH_VARARGS, ""},
  {NULL, NULL, 0, NULL}
};

// init(modulename)
PyMODINIT_FUNC initlib_LBScan_c(void)
{
  (void)Py_InitModule("lib_LBScan_c", methods);  // "modulename"
  import_array();
}
