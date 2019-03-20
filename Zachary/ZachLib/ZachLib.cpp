/*

    Implements some performance critical portions of the Zachary package in PyVPT

*/

#include "ZachLib.h"
#include "PyExtLib.h"

/******************************** BEGIN MAIN PACKAGE BODY ************************************/

/*
 * ZachLib_UnevenFiniteDifferenceWeights
 *
 * we use the formula from
       https://pdfs.semanticscholar.org/8bf5/912bde884f6bd4cfb4991ba3d077cace94c0.pdf
 * and fill a preallocated numpy array
 */

int ind2d(int i, int j, int n, int m) {
    return m*i + j;
}
FUNCWITHARGS(ZachLib_UnevenFiniteDifferenceWeights) {

    Py_ssize_t m;
    double z;
    PyObject *grid;
    PARSEARGS("ndO", &m, &z, &grid);

    PyObject *xArray = _ArrayAsType(grid, "float64");
    CHECKNULL(xArray);
    double *x = _GetDataArray<double>(xArray);
    CHECKCLEAN(x, xArray);

    // build the appropriate nxm array to fill
    Py_ssize_t n = PyObject_Length(xArray) - 1;
    int dims[2] = {n+1, m+1};
    PyObject *cArray = _CreateArray(2, dims, "zeros");
    double *c = _GetDataArray<double>(cArray);
//    PyObject *strides = PyObject_GetAttrString(cArray, "strides");
//    _DebugPrintObject(0, strides);
//    CLEANUP(strides);
    // basic implementation just copied from the original source
    double dxProdOld, dxProd, dx, dz, dzOld;
    dxProdOld = 1.;
    dz = x[0] - z;
    c[0] = 1.;
    double weight = 0.;
    int mn;
    int ind1, ind2, ind3;
    int nrow = n+1;
    int ncols = m+1;
    for ( int i = 1; i<=n; i++) {// do i=1,n
        if (i < m) { mn = i; } else { mn = m; }; // mn = min(i, m)
        dxProd = 1.; // c2 = 1.0d0
        dzOld = dz; // c5 = c4
        dz = x[i] - z; // c4 = x(i)-z
        for ( int j = 0; j<i; j++ ) { // do 40 j=0,i-1
            dx = x[i] - x[j]; // c3 = x(i)-x(j)
            dxProd = dxProd * dx; // c2 = c2*c3
            if (j == i-1) { // if (j.eq.i-1) then
                for (int k = mn; k>0; k--){ // do 20 k=mn,1,-1
                    ind1 = ind2d(i,   k,   nrow, ncols); // (i,k)
                    ind2 = ind2d(i-1, k-1, nrow, ncols); // (i-1,k-1)
                    ind3 = ind2d(i-1, k,   nrow, ncols); // (i-1,k)
                    weight = (k*c[ind2] - dzOld*c[ind3])*dxProdOld/dxProd;
                    c[ind1] = weight; // c(i,k) = c1*(k*c(i-1,k-1)-c5*c(i-1,k))/c2
                }
                ind1 = ind2d(i,   0, nrow, ncols); // (i,0)
                ind2 = ind2d(i-1, 0, nrow, ncols); // (i-1, 0)
                c[ind1] = -dzOld*c[ind2]*dxProdOld/dxProd; // c(i,0) = -c1*c5*c(i-1,0)/c2
            }

            for (int k = mn; k>0; k--){ // do 30 k=mn,1,-1
                ind1 = ind2d(j, k,   nrow, ncols); // (j, k)
                ind2 = ind2d(j, k-1, nrow, ncols); // (j, k-1)
                weight = (dz*c[ind1]-k*c[ind2])/dx;
                c[ind1] = weight; // c(j,k) = (c4*c(j,k)-k*c(j,k-1))/c3
            }
            ind1 = ind2d(j, 0, nrow, ncols); // (j,0)
            weight = dz*c[ind1]/dx;
            c[ind1] = weight;  //  c(j,0) = c4*c(j,0)/c3
        }

        dxProdOld = dxProd; //  c1 = c2

    }

    // return final numpy array
    return cArray;

}

/******************************** MODULE INIT BITS ********************************/

static PyMethodDef ZachLibMethods[] = {
    {"UnevenFiniteDifferenceWeights", ZachLib_UnevenFiniteDifferenceWeights, METH_VARARGS, ""},
};

const char ZachLib_doc[] = "ZachLib is a layer for doing performance critical stuff";

static struct PyModuleDef ZachLibModule = {
    PyModuleDef_HEAD_INIT,
    "ZachLib",   /* name of module */
    ZachLib_doc, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    ZachLibMethods
};

PyMODINIT_FUNC PyInit_ZachLib(void)
{
    return PyModule_Create(&ZachLibModule);
}


