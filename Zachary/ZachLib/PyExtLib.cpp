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

PyObject *_CreateArray(int depth, int *dims, const char *ctor) {
    PyObject *array_module = PyImport_ImportModule("numpy");
    CHECKNULL(array_module);
    PyObject *builder = PyObject_GetAttrString(array_module, ctor);
    CHECKCLEAN(builder, array_module);
    PyObject *dimObj = PyList_New(depth);
    for (int j = 0; j<depth; j++){
        PyList_SetItem(dimObj, j, Py_BuildValue("i", dims[j]));
    }
    PyObject *cArray = PyObject_CallFunction(builder, "O", dimObj);
    CHECKCLEAN(cArray, builder, array_module);
    CLEANUP(array_module, builder);

    return cArray;
}


