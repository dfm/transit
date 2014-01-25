#include <Python.h>
#include <numpy/arrayobject.h>

#include "driver.h"

/* Docstrings */
static char doc[] = "Compute transit models.\n";

static PyObject *transit_ldlc(PyObject *self, PyObject *args)
{
    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef transit_methods[] = {
    {"ldlc", transit_ldlc, METH_VARARGS, ""},
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
