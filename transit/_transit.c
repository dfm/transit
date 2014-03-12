#include <Python.h>
#include <numpy/arrayobject.h>

#include "driver.h"

#define PARSE_ARRAY(o) (PyArrayObject*) PyArray_FROM_OTF(o, NPY_DOUBLE, NPY_IN_ARRAY)

/* Docstrings */
static char doc[] = "Compute transit models.\n";

static PyObject *transit_ldlc_kepler(PyObject *self, PyObject *args)
{
    int maxdepth;
    double mstar, rstar, mu1, mu2, texp, tol;
    PyObject *t_obj, *occ_obj, *mp_obj, *r_obj, *a_obj, *t0_obj, *e_obj,
             *pomega_obj, *ix_obj, *iy_obj;

    // Parse the input arguments.
    if (!PyArg_ParseTuple(args, "OddddOOOOOOOOOddi",
                          &t_obj, &mu1, &mu2, &mstar, &rstar, &occ_obj,
                          &mp_obj, &r_obj, &a_obj, &t0_obj, &e_obj,
                          &pomega_obj, &ix_obj, &iy_obj,
                          &texp, &tol, &maxdepth))
        return NULL;

    // Decode the numpy arrays.
    PyArrayObject *t_array = PARSE_ARRAY(t_obj),
                  *occ_array = PARSE_ARRAY(occ_obj),
                  *mp_array = PARSE_ARRAY(mp_obj),
                  *r_array = PARSE_ARRAY(r_obj),
                  *a_array = PARSE_ARRAY(a_obj),
                  *t0_array = PARSE_ARRAY(t0_obj),
                  *e_array = PARSE_ARRAY(e_obj),
                  *pomega_array = PARSE_ARRAY(pomega_obj),
                  *ix_array = PARSE_ARRAY(ix_obj),
                  *iy_array = PARSE_ARRAY(iy_obj);
    if (t_array == NULL || occ_array == NULL || mp_array == NULL ||
            r_array == NULL || a_array == NULL || t0_array == NULL ||
            e_array == NULL || pomega_array == NULL || ix_array == NULL ||
            iy_array == NULL)
        goto fail;

    // Check the dimensions.
    int n = (int) PyArray_DIM(t_array, 0),
        np = (int) PyArray_DIM(mp_array, 0);
    if (
            (int)PyArray_DIM(occ_array, 0) != np    ||
            (int)PyArray_DIM(r_array, 0) != np      ||
            (int)PyArray_DIM(a_array, 0) != np      ||
            (int)PyArray_DIM(t0_array, 0) != np     ||
            (int)PyArray_DIM(e_array, 0) != np      ||
            (int)PyArray_DIM(pomega_array, 0) != np ||
            (int)PyArray_DIM(ix_array, 0) != np     ||
            (int)PyArray_DIM(iy_array, 0) != np
    ) {
        PyErr_SetString(PyExc_ValueError, "Dimension mismatch");
        goto fail;
    }

    // Get pointers to the input data.
    double *t = PyArray_DATA(t_array),
           *occ = PyArray_DATA(occ_array),
           *mp = PyArray_DATA(mp_array),
           *r = PyArray_DATA(r_array),
           *a = PyArray_DATA(a_array),
           *t0 = PyArray_DATA(t0_array),
           *e = PyArray_DATA(e_array),
           *pomega = PyArray_DATA(pomega_array),
           *ix = PyArray_DATA(ix_array),
           *iy = PyArray_DATA(iy_array);

    // Allocate the flux array.
    npy_intp dim[1] = {n};
    PyArrayObject *lam_array = (PyArrayObject*)PyArray_SimpleNew(1, dim, NPY_DOUBLE);
    if (lam_array == NULL) {
        Py_XDECREF(lam_array);
        goto fail;
    }

    // Compute the model light curve.
    double *lam = PyArray_DATA(lam_array);
    int info = ldlc_kepler(mu1, mu2, mstar, rstar, np, occ, mp, r, a, t0, e, pomega, ix, iy, texp, tol, maxdepth, n, t, lam);

    // Clean up a bit.
    Py_DECREF(t_array);
    Py_DECREF(occ_array);
    Py_DECREF(mp_array);
    Py_DECREF(r_array);
    Py_DECREF(a_array);
    Py_DECREF(t0_array);
    Py_DECREF(e_array);
    Py_DECREF(pomega_array);
    Py_DECREF(ix_array);
    Py_DECREF(iy_array);

    if (info) {
        Py_DECREF(lam_array);
        PyErr_SetString(PyExc_RuntimeError, "Orbit solve failed.");
        return NULL;
    }

    return (PyObject*)lam_array;

fail:

    Py_XDECREF(t_array);
    Py_XDECREF(occ_array);
    Py_XDECREF(mp_array);
    Py_XDECREF(r_array);
    Py_XDECREF(a_array);
    Py_XDECREF(t0_array);
    Py_XDECREF(e_array);
    Py_XDECREF(pomega_array);
    Py_XDECREF(ix_array);
    Py_XDECREF(iy_array);
    return NULL;
}

static PyObject *transit_ldlc_simple(PyObject *self, PyObject *args)
{
    int maxdepth;
    double mu1, mu2, p, t0, tau, ror, b, texp, tol;
    PyObject *t_obj;

    // Parse the input arguments.
    if (!PyArg_ParseTuple(args, "Odddddddddi",
                          &t_obj, &mu1, &mu2, &p, &t0, &tau, &ror, &b,
                          &texp, &tol, &maxdepth))
        return NULL;

    // Decode the numpy arrays.
    PyArrayObject *t_array = PARSE_ARRAY(t_obj);
    if (t_array == NULL)
        return NULL;

    int n = (int) PyArray_DIM(t_array, 0);
    double *t = PyArray_DATA(t_array);

    // Allocate the flux array.
    npy_intp dim[1] = {n};
    PyArrayObject *lam_array = (PyArrayObject*)PyArray_SimpleNew(1, dim, NPY_DOUBLE);
    if (lam_array == NULL) {
        Py_DECREF(t_array);
        Py_XDECREF(lam_array);
        return NULL;
    }

    // Compute the model light curve.
    double *lam = PyArray_DATA(lam_array);
    int info = ldlc_simple(mu1, mu2, p, t0, tau, ror, b, texp, tol, maxdepth, n, t, lam);

    // Clean up a bit.
    Py_DECREF(t_array);

    if (info) {
        Py_DECREF(lam_array);
        PyErr_SetString(PyExc_RuntimeError, "Orbit solve failed.");
        return NULL;
    }

    return (PyObject*)lam_array;
}

static PyMethodDef transit_methods[] = {
    {"ldlc_kepler", transit_ldlc_kepler, METH_VARARGS, ""},
    {"ldlc_simple", transit_ldlc_simple, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3

#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

struct module_state {
    PyObject *error;
};

static int transit_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int transit_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_transit",
    doc,
    sizeof(struct module_state),
    transit_methods,
    NULL,
    transit_traverse,
    transit_clear,
    NULL
};

#define INITERROR return NULL

PyObject *PyInit__transit(void)
#else

#define INITERROR return

void init_transit(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule3("_transit", transit_methods, doc);
#endif

    if (module == NULL)
        INITERROR;
    import_array();

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
