#include "Python.h"
#include "arrayobject.h"

static PyObject *_getMultipleDists( PyObject * self, PyObject * args )
{
	PyArrayObject *ary_QList, *ary_QConf;
	if ( ! PyArg_ParseTuple( args,"iiiOOOf", &ary_QList, &ary_QConf ) )
	{
		return NULL;
	}

}
