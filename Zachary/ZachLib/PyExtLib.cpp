/*
    Implements some common useful utilities in writing Python extensions
*/

#include "PyExtLib.h"

/******************************** ARRAY MANIPULATION HELPERS ********************************/

Py_buffer _GetDataBuffer(PyObject *data) {

    Py_buffer view;
    PyObject_GetBuffer(data, &view, PyBUF_CONTIG_RO);
    return view;
}

//template<typename T>
//T *_GetDataBufferArray(Py_buffer *view) {

//    T *c_data;
//    if ( view == NULL ) return NULL;
//    c_data = (T *) view->buf;
//    if (c_data == NULL) {
//        PyBuffer_Release(view);
//    }
//    return c_data;

//}

//template<typename T>
//T *_GetDataArray(PyObject *data) {
//    Py_buffer view = _GetDataBuffer(data);
//    T *array = _GetDataBufferArray<T>(&view);
//    CHECKNULL(array);
//    return array;
//}

void _CopyDataBuffer(PyObject *data, void *buff, long len, int dsize) {

    Py_buffer view_obj = _GetDataBuffer(data);
    Py_buffer *view = &view_obj;
    if (len) memcpy(view->buf, buff, len*dsize);

}

//template<typename T>
//void _CopyDataBuffer(PyObject *data, void *buff, long len) {
//    _CopyDataBuffer(data, buff, len, sizeof(T));
//}

//template<typename T>
//void _SetDataBuffer(PyObject *data, void *buff) {

//    Py_buffer view_obj = _GetDataBuffer(data);
//    Py_buffer *view = &view_obj;
//    view->buf = (T *) buff;

//}

int _DebugPrint(int level, const char *fmt, ...) {
    if (level <= _DEBUG_LEVEL) {
        va_list args;
        va_start(args, fmt);

        int res = vprintf(fmt, args);
        printf("\n");

        fflush(stdout);

        return res;
    } else {
        return 0;
    }
}
int _DebugMessage(int level, const char *msg) {
    return _DebugPrint(level, "%s", msg);
}
const char *_Repr(PyObject *o, PyObject *repr) {
    PyObject *tmp = PyObject_Repr(o);
    PyObject *enc = PyUnicode_AsEncodedString(tmp, "utf-8", "strict");
    if ( enc == NULL) {
        Py_XDECREF(tmp);
        return NULL;
    }
    const char *str =  PyBytes_AsString(enc);
    Py_XDECREF(enc);
    Py_XDECREF(tmp);

    return str;
}
int _DebugPrintObject(int lvl, PyObject *o){
    PyObject *repr = NULL;
    const char * buff=_Repr(o, repr);
    int res = _DebugPrint(lvl, buff);
    Py_XDECREF(repr);
    return res;
}

PyObject *_ArrayAsType(PyObject *array, const char *type) {

    PyObject *array_module = PyImport_ImportModule("numpy");
    CHECKNULL(array_module)
    PyObject *astype = PyObject_GetAttrString(array, "astype");
    CHECKCLEAN(astype, array_module)
    PyObject *targetType = PyObject_GetAttrString(array_module, type);
    CHECKCLEAN(targetType, astype, array_module)
    PyObject *xArray = PyObject_CallFunction(astype, "O", targetType);
    CHECKCLEAN(xArray, targetType, astype, array_module)

    return xArray;

}

PyObject *_CreateFromNumPy(const char *ctor, PyObject *args, PyObject *kwargs) {

    PyObject *array_module = PyImport_ImportModule("numpy");
    CHECKNULL(array_module);
    PyObject *builder = PyObject_GetAttrString(array_module, ctor);
    CHECKCLEAN(builder, array_module);
    PyObject *cArray = PyObject_Call(builder, args, kwargs);
    CHECKCLEAN(cArray, builder, array_module);
    CLEANUP(builder, array_module);

    return cArray;
}
PyObject *_CreateFromNumPy(const char *ctor, PyObject *args, const char *dtype) {

    PyObject *array_module = PyImport_ImportModule("numpy");
    CHECKNULL(array_module);
    PyObject *typeObj = PyObject_GetAttrString(array_module, dtype);
    CHECKCLEAN(typeObj, array_module);
    PyObject *kwargs = Py_BuildValue("{s:O}", "dtype", typeObj);
    CHECKCLEAN(kwargs, typeObj, array_module);

    PyObject *cArray = _CreateFromNumPy(ctor, args, kwargs);
    CHECKCLEAN(cArray, kwargs, typeObj, array_module);
    CLEANUP(kwargs, typeObj, array_module);

    return cArray;
}

PyObject *_CreateIdentity(int m, const char *dtype) {

    PyObject *args = Py_BuildValue("(i)", m);
    CHECKNULL(args);
    PyObject *cArray = _CreateFromNumPy("eye", args, dtype);
    CHECKCLEAN(cArray, args);
    CLEANUP(args);

    return cArray;
}
PyObject *_CreateIdentity(int m) {
    return _CreateIdentity(m, "int8");
}
PyObject *_CreateArray(int depth, int *dims, const char *ctor, const char *dtype) {

    PyObject *dimObj = PyList_New(depth);
    CHECKNULL(dimObj);
    for (int j = 0; j<depth; j++){
        PyList_SetItem(dimObj, j, Py_BuildValue("i", dims[j]));
    }
    PyObject *args = Py_BuildValue("(O)", dimObj);
    PyObject *cArray = _CreateFromNumPy(ctor, args, dtype);
    CHECKCLEAN(cArray, args, dimObj);
    CLEANUP(args, dimObj);

    return cArray;
}
PyObject *_CreateArray(int depth, int *dims, const char *ctor) {
    return _CreateArray(depth, dims, ctor, "float64");
}


