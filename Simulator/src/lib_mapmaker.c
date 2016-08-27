#include "/usr/include/python2.7/Python.h"
#include "/usr/lib64/python2.7/site-packages/numpy/core/include/numpy/arrayobject.h"
#include "/usr/lib64/python2.7/site-packages/numpy/core/include/numpy/arrayscalars.h"
#include <bits/nan.h> 

// #include "/usr/common/usg/python/2.7.1/include/python2.7/Python.h"
// #include "/usr/common/usg/python/numpy/1.6.1/lib/python/numpy/core/include/numpy/arrayobject.h"
// #include "/usr/common/usg/python/numpy/1.6.1/lib/python/numpy/core/include/numpy/arrayscalars.h"

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

int not_doublevector(PyArrayObject *vec){
  if (vec->descr->type_num != NPY_DOUBLE || vec->nd != 1)  {
    PyErr_SetString(PyExc_ValueError,
		    "In not_doublevector: array must be of type Float and 1 dimensional (n).");
    return 1;  }
  return 0;
}

//##############################################################################################

static PyObject *MapProjection_forLoopInC_part(PyObject *self, PyObject *args)
{
  //  printf("[MapProjection_forLoopInC_part in lib_mapmaker.c] START ");

  int nb;
  PyArrayObject *i_hpx, *sum_dat, *dif_dat, *alpha, *beta, *p_weight, *n_weight, *ind_all;
  int dims[2];
  double *c_i_hpx, *c_sum_dat, *c_dif_dat, *c_alpha, *c_beta, *c_p_weight, *c_n_weight, *c_ind_all;
  PyArrayObject *mapout;
  double **c_mapout;


  /* Parse tuples separately since args will differ between C fcns */
  if (!PyArg_ParseTuple(args, "iO!O!O!O!O!O!O!O!", 
			&nb,
			&PyArray_Type, &i_hpx,
			&PyArray_Type, &sum_dat,
			&PyArray_Type, &dif_dat,
			&PyArray_Type, &alpha,
			&PyArray_Type, &beta,
			&PyArray_Type, &p_weight,
			&PyArray_Type, &n_weight,
			&PyArray_Type, &ind_all))  return NULL;
  if (NULL == i_hpx)  return NULL;
  if (NULL == sum_dat)  return NULL;
  if (NULL == dif_dat)  return NULL;
  if (NULL == alpha)  return NULL;
  if (NULL == beta)  return NULL;
  if (NULL == p_weight)  return NULL;
  if (NULL == n_weight)  return NULL;
  if (NULL == ind_all)  return NULL;
  
  /* Check that objects are 'double' type and vectors
     Not needed if python wrapper function checks before call to this routine */
  if (not_doublevector(i_hpx)) return NULL;
  if (not_doublevector(sum_dat)) return NULL;
  if (not_doublevector(dif_dat)) return NULL;
  if (not_doublevector(alpha)) return NULL;
  if (not_doublevector(beta)) return NULL;
  if (not_doublevector(p_weight)) return NULL;
  if (not_doublevector(n_weight)) return NULL;
  if (not_doublevector(ind_all)) return NULL;
  
  /* Get vector dimension. */
  //  printf("%d \n", nb);
  dims[0] = 10;
  dims[1] = nb;
  mapout = (PyArrayObject *) PyArray_FromDims(2,dims,NPY_DOUBLE);
  
  /* Change contiguous arrays into C * arrays   */
  c_i_hpx = pyvector_to_Carrayptrs(i_hpx);
  c_sum_dat = pyvector_to_Carrayptrs(sum_dat);
  c_dif_dat = pyvector_to_Carrayptrs(dif_dat);
  c_alpha = pyvector_to_Carrayptrs(alpha);
  c_beta = pyvector_to_Carrayptrs(beta);
  c_p_weight = pyvector_to_Carrayptrs(p_weight);
  c_n_weight = pyvector_to_Carrayptrs(n_weight);
  c_ind_all = pyvector_to_Carrayptrs(ind_all);
  
  c_mapout = pymatrix_to_Carrayptrs(mapout);
  
  //  printf("%d %d %d %d %d \n", nb, i_hpx->dimensions[0], i_hpx->dimensions[1], ind_all->dimensions[0], ind_all->dimensions[1]);
  int i, ii, nb_tod=i_hpx->dimensions[0];
  //  printf("%d %d\n", nb, nb_tod);
  //  printf("[MapProjection_forLoopInC_part in lib_mapmaker.c] %d %d %d %d %d %d\n", 
  //	 nb, i_hpx->dimensions[0], i_hpx->dimensions[1], ind_all->dimensions[0], ind_all->dimensions[1], nb_tod);
  for (i=0;i<nb_tod;i++){
    //    printf("%f %d \n", c_ind_all[(int)c_i_hpx[i]], (int)c_ind_all[(int)c_i_hpx[i]] );
    ii = (int)c_ind_all[(int)c_i_hpx[i]];
    //    printf("%d %d %d \n", i, (int)c_i_hpx[i], ii);
    c_mapout[0][ii] += c_p_weight[i]*c_sum_dat[i];
    c_mapout[1][ii] += c_p_weight[i];
    c_mapout[2][ii] += c_alpha[i] * c_n_weight[i];
    c_mapout[3][ii] +=  c_beta[i] * c_n_weight[i];
    c_mapout[4][ii] += c_alpha[i] *   c_alpha[i] * c_n_weight[i];
    c_mapout[5][ii] +=  c_beta[i] *    c_beta[i] * c_n_weight[i];
    c_mapout[6][ii] += c_alpha[i] *    c_beta[i] * c_n_weight[i];
    c_mapout[7][ii] += c_alpha[i] * c_dif_dat[i] * c_n_weight[i];
    c_mapout[8][ii] +=  c_beta[i] * c_dif_dat[i] * c_n_weight[i];
    c_mapout[9][ii] += 1.;
  }

  //  free_Carrayptrs(c_i_hpx);
  //  free_Carrayptrs(c_sum_dat);
  //  free_Carrayptrs(c_dif_dat);
  //  free_Carrayptrs(c_alpha);
  // free_Carrayptrs(c_beta);
  //  free_Carrayptrs(c_p_weight);
  //  free_Carrayptrs(c_n_weight);
  //  free_Carrayptrs(c_ind_all);
  
  //  printf("[MapProjection_forLoopInC_part in lib_mapmaker.c] END \n");
  //  return Py_BuildValue("i", 1);
  return PyArray_Return(mapout);
}

//##############################################################################################

static PyObject *MapProjection_forLoopInC(PyObject *self, PyObject *args)
{
  //  printf("[MapProjection_forLoopInC_part in lib_mapmaker.c] START ");

  int nb;
  PyArrayObject *i_hpx, *sum_dat, *dif_dat, *alpha, *beta, *p_weight, *n_weight;
  int dims[2];
  double *c_i_hpx, *c_sum_dat, *c_dif_dat, *c_alpha, *c_beta, *c_p_weight, *c_n_weight;
  PyArrayObject *mapout;
  double **c_mapout;


  /* Parse tuples separately since args will differ between C fcns */
  if (!PyArg_ParseTuple(args, "iO!O!O!O!O!O!O!", 
			&nb,
			&PyArray_Type, &i_hpx,
			&PyArray_Type, &sum_dat,
			&PyArray_Type, &dif_dat,
			&PyArray_Type, &alpha,
			&PyArray_Type, &beta,
			&PyArray_Type, &p_weight,
			&PyArray_Type, &n_weight))  return NULL;
  if (NULL == i_hpx)  return NULL;
  if (NULL == sum_dat)  return NULL;
  if (NULL == dif_dat)  return NULL;
  if (NULL == alpha)  return NULL;
  if (NULL == beta)  return NULL;
  if (NULL == p_weight)  return NULL;
  if (NULL == n_weight)  return NULL;
  
  /* Check that objects are 'double' type and vectors
     Not needed if python wrapper function checks before call to this routine */
  if (not_doublevector(i_hpx)) return NULL;
  if (not_doublevector(sum_dat)) return NULL;
  if (not_doublevector(dif_dat)) return NULL;
  if (not_doublevector(alpha)) return NULL;
  if (not_doublevector(beta)) return NULL;
  if (not_doublevector(p_weight)) return NULL;
  if (not_doublevector(n_weight)) return NULL;
  
  /* Get vector dimension. */
  //  printf("%d \n", nb);
  dims[0] = 11;
  dims[1] = nb;
  mapout = (PyArrayObject *) PyArray_FromDims(2,dims,NPY_DOUBLE);
  
  /* Change contiguous arrays into C * arrays   */
  c_i_hpx = pyvector_to_Carrayptrs(i_hpx);
  c_sum_dat = pyvector_to_Carrayptrs(sum_dat);
  c_dif_dat = pyvector_to_Carrayptrs(dif_dat);
  c_alpha = pyvector_to_Carrayptrs(alpha);
  c_beta = pyvector_to_Carrayptrs(beta);
  c_p_weight = pyvector_to_Carrayptrs(p_weight);
  c_n_weight = pyvector_to_Carrayptrs(n_weight);
  
  c_mapout = pymatrix_to_Carrayptrs(mapout);
  
  //  printf("%d %d %d %d %d \n", nb, i_hpx->dimensions[0], i_hpx->dimensions[1], ind_all->dimensions[0], ind_all->dimensions[1]);
  int i, ii, nb_tod=i_hpx->dimensions[0];
  //  printf("%d %d\n", nb, nb_tod);
  //  printf("[MapProjection_forLoopInC_part in lib_mapmaker.c] %d %d %d %d %d %d\n", 
  //	 nb, i_hpx->dimensions[0], i_hpx->dimensions[1], ind_all->dimensions[0], ind_all->dimensions[1], nb_tod);
  for (i=0;i<nb_tod;i++){
    //    printf("%f %d \n", c_ind_all[(int)c_i_hpx[i]], (int)c_ind_all[(int)c_i_hpx[i]] );
    ii = (int)c_i_hpx[i];
    //    printf("%d %d %d \n", i, (int)c_i_hpx[i], ii);
    //    c_mapout[0][ii] += c_p_weight[i]*c_sum_dat[i];
    //    c_mapout[1][ii] += c_p_weight[i];
    c_mapout[0][ii] += c_sum_dat[i];  // In
    c_mapout[1][ii] += 1.;  // Id
    c_mapout[2][ii] += c_alpha[i]; // * c_n_weight[i];, A
    c_mapout[3][ii] +=  c_beta[i]; //  * c_n_weight[i];, B
    c_mapout[4][ii] += c_alpha[i] *   c_alpha[i]; // * c_n_weight[i];, AA
    c_mapout[5][ii] +=  c_beta[i] *    c_beta[i]; // * c_n_weight[i];, BB
    c_mapout[6][ii] += c_alpha[i] *    c_beta[i]; // * c_n_weight[i];, AB
    c_mapout[7][ii] += c_alpha[i] * c_dif_dat[i]; // * c_n_weight[i];, Ad
    c_mapout[8][ii] +=  c_beta[i] * c_dif_dat[i]; // * c_n_weight[i];, Bd
    c_mapout[9][ii] += c_dif_dat[i];
    c_mapout[10][ii] += 1.; // H
  }
  
  //  free_Carrayptrs(c_i_hpx);
  //  free_Carrayptrs(c_sum_dat);
  //  free_Carrayptrs(c_dif_dat);
  //  free_Carrayptrs(c_alpha);
  // free_Carrayptrs(c_beta);
  //  free_Carrayptrs(c_p_weight);
  //  free_Carrayptrs(c_n_weight);
  //  free_Carrayptrs(c_ind_all);
  
  //  printf("[MapProjection_forLoopInC_part in lib_mapmaker.c] END \n");
  //  return Py_BuildValue("i", 1);
  return PyArray_Return(mapout);
}

//##############################################################################################

static PyObject *boreptg2pixptg_c(PyObject *self, PyObject *args)
{
  int i;

  /* Parse tuples separately since args will differ between C fcns */
  int nb;
  double fp_x, fp_y, fp_polang;
  PyArrayObject *b_ra, *b_dec, *b_dk, *b_hwp;
  double *c_b_ra, *c_b_dec, *c_b_dk, *c_b_hwp;
  double radeg = (180./M_PI);

  double beam_amp = (2./60./180.*M_PI);
  double theta_rad, phi_rad;
  double d_x, d_y, d_z, x_bt, y_bt, z_bt, x_bp, y_bp, z_bp;
  double x, y, z, x_tmp;
  double polang_beta, chi;

  if (!PyArg_ParseTuple(args, "iO!O!O!O!ddd", 
			&nb,
			&PyArray_Type, &b_ra,
			&PyArray_Type, &b_dec,
			&PyArray_Type, &b_dk,
			&PyArray_Type, &b_hwp,
			&fp_x,
			&fp_y,
			&fp_polang))  return NULL;
  if (NULL == b_ra)  return NULL;
  if (NULL == b_dec)  return NULL;
  if (NULL == b_dk)  return NULL;
  if (NULL == b_hwp)  return NULL;
  
  /* Check that objects are 'double' type and vectors
     Not needed if python wrapper function checks before call to this routine */
  if (not_doublevector(b_ra)) return NULL;
  if (not_doublevector(b_dec)) return NULL;
  if (not_doublevector(b_dk)) return NULL;
  if (not_doublevector(b_hwp)) return NULL;

  int dims[3];
  dims[0] = 3;
  dims[1] = nb;
  PyArrayObject *pixang_out;
  double **c_pixang_out;
  pixang_out = (PyArrayObject *) PyArray_FromDims(2,dims,NPY_DOUBLE);
  c_pixang_out = pymatrix_to_Carrayptrs(pixang_out);

  /* Change contiguous arrays into C * arrays   */
  c_b_ra = pyvector_to_Carrayptrs(b_ra);
  c_b_dec = pyvector_to_Carrayptrs(b_dec);
  c_b_dk = pyvector_to_Carrayptrs(b_dk);
  c_b_hwp = pyvector_to_Carrayptrs(b_hwp);
 
  double theta_b, phi_b, dk, theta_fp, phi_fp;
  double theta_p, phi_p;
  double x1_, y1_, z1_, x2_, y2_, z2_;
  double vec_theta_x, vec_theta_y, vec_theta_z;
  double vec_phi_x, vec_phi_y, vec_phi_z;
  double beta;
  double tmp, tmp_prev=0.;

  for (i=0;i<nb;i++){
    //    printf("fp_polang %lf \n", fp_polang);
    theta_b =  M_PI/2.-c_b_dec[i];
    phi_b = c_b_ra[i];
    dk = c_b_dk[i];
    theta_fp = sqrt(pow(fp_x,2)+pow(fp_y,2));
    phi_fp = atan2(fp_y,fp_x);

    // compute the pixel pointing    
    //    if (i==0){
    //      printf("%lf %lf %lf %lf %lf \n", theta_b*radeg, theta_fp*radeg, dk*radeg, phi_fp*radeg, phi_b*radeg);
    //    }
    theta_p = acos( -sin(theta_b)*sin(theta_fp)*cos(dk+phi_fp)+cos(theta_b)*cos(theta_fp) );
    phi_p = phi_b + atan2( (sin(theta_fp)*sin(dk+phi_fp)), (cos(theta_b)*sin(theta_fp)*cos(dk+phi_fp)+sin(theta_b)*cos(theta_fp)) );
    //    if (i==0){
    //      printf("%lf %lf \n", theta_p*radeg, phi_p*radeg);
    //    }

    // polarization angle at the focal plane coordinate
    x1_ = cos(fp_polang)*cos(theta_fp)*cos(phi_fp+dk)-sin(fp_polang)*sin(phi_fp+dk);
    y1_ = cos(fp_polang)*cos(theta_fp)*sin(phi_fp+dk)+sin(fp_polang)*cos(phi_fp+dk);
    z1_ = -cos(fp_polang)*sin(theta_fp);

    // rotation for the polarization angle
    x2_ =  cos(theta_b)*x1_+sin(theta_b)*z1_;
    y2_ = y1_;
    z2_ = -sin(theta_b)*x1_+cos(theta_b)*z1_;

    x1_ = cos(phi_b)*x2_-sin(phi_b)*y2_;
    y1_ = sin(phi_b)*x2_+cos(phi_b)*y2_;
    z1_ = z2_;

    // compute d_theta, d_phi at the pixel pointing
    vec_theta_x = cos(theta_p)*cos(phi_p);
    vec_theta_y = cos(theta_p)*sin(phi_p);
    vec_theta_z = -sin(theta_p);

    vec_phi_x = -sin(phi_p);
    vec_phi_y = cos(phi_p);
    vec_phi_z = 0.;

    // compute the rotated polarization angle w.r.t. d_theta and d_phi
    // and compute the polarization angle with healpix convention
    tmp = acos( x1_*vec_theta_x + y1_*vec_theta_y + z1_*vec_theta_z);
    if (isnan(tmp)) tmp = tmp_prev;

    c_pixang_out[2][i] = tmp;
    beta = acos(  x1_*vec_phi_x + y1_*vec_phi_y + z1_*vec_phi_z);
    if (beta>90./radeg){
      c_pixang_out[2][i] = - c_pixang_out[2][i];
    }
    c_pixang_out[0][i] = phi_p;
    c_pixang_out[1][i] = M_PI/2.-theta_p;
    tmp_prev = tmp;
    
  }
  
  return PyArray_Return(pixang_out);

}

//##############################################################################################

static PyObject *TestInC(PyObject *self, PyObject *args)
{
  //  printf("[MapProjection_forLoopInC_part in lib_mapmaker.c] START ");
  int i;
  PyArrayObject *array;
  double *c_array;

  /* Parse tuples separately since args will differ between C fcns */
  if (!PyArg_ParseTuple(args, "O!", 
			&PyArray_Type, &array))  return NULL;
  if (NULL == array)  return NULL;
  if (not_doublevector(array)) return NULL;
  c_array = pyvector_to_Carrayptrs(array);
  for (i=0;i<10;i++) {
    printf("array in C %lf \n\n",c_array[i]);
  }
  return PyArray_Return(array);
}

//##############################################################################################

static PyMethodDef methods[] = {
  {"MapProjection_forLoopInC_part", MapProjection_forLoopInC_part, METH_VARARGS, ""},
  {"MapProjection_forLoopInC", MapProjection_forLoopInC, METH_VARARGS, ""},
  {"boreptg2pixptg_c", boreptg2pixptg_c, METH_VARARGS, ""},
  {"TestInC", TestInC, METH_VARARGS, ""},
  {NULL, NULL, 0, NULL}
};

// init(modulename)
PyMODINIT_FUNC init_lib_mapmaker(void)
{
  (void)Py_InitModule("_lib_mapmaker", methods);  // "modulename"
  import_array();
}
