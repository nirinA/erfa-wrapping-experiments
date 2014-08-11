'''erfamodule generator
'''
from __future__ import with_statement

try:
    next
except NameError:
    def next(n):
        return n.next()

class ErfaWrapp:
    def __init__(self):
        '''read erfa.c and create a dict with the name
of each function, its args, type and doc'''
        self.functions = {}
        with open('erfa.c') as f:
            g = (l for l in  f.readlines())
            for l in g:
                if l.startswith(('void', 'int', 'double')):
                    statement = l
                    while not statement.endswith(')\n'):
                        statement += next(g)
                    func, sep, args = statement.partition('(')
                    type_, name = func.split()
                    args = args.rstrip(')\n').replace('\n', ' ').split(',')
                    doc = ''
                    while not doc.endswith('*/\n'):
                        doc += next(g)
                    self.functions.update({name:{'type':type_,
                                                 'args':args,
                                                 'doc':doc}})
        for n in self.functions.keys():
            self.parse_doc(n)

    def parse_doc(self, name):
        '''update function dict with lists of input, output arguments and returned value'''
        given = []
        returned = []
        returned_value = []
        g = (l for l in self.functions[name]['doc'].split('\n'))
        for l in g:
            if 'Given' in l:
                given.append(next(g))
                p = next(g)
                while p != '**':
                    given.append(p)
                    p = next(g)
            elif 'Returned' in l:
                if '(function value):' in l:
                    returned_value.append(next(g))
                    p = next(g)
                    while p != '**':
                        returned_value.append(p)
                        p = next(g)
                else:
                    returned.append(next(g))
                    p = next(g)
                    while p != '**':
                        returned.append(p)
                        p = next(g)
        self.functions[name]['returned_value'] = returned_value
        args_in = []
        for i in given:
            n = i.split()[1]
            if ',' in n:
                args_in.extend(n.split(','))
            else:
                args_in.append(n)
        args_out = []
        for i in returned:
            n = i.split()[1]
            if ',' in n:
                args_out.extend(n.split(','))
            else:
                args_out.append(n)        
        self.functions[name]['args_in'] = args_in
        self.functions[name]['args_out'] = args_out

## Templates
header_templates = '''#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <math.h>
#include "numpy/npy_3kcompat.h"
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"
#include "numpy/npy_math.h"
#include "erfa.h"

static PyObject *_erfaError;

static PyObject *
_erfa_bi00(PyObject *self)
{
    double dpsibi, depsbi, dra;
    eraBi00(&dpsibi, &depsbi, &dra);
    if (!PyErr_Occurred()) {
        return Py_BuildValue("ddd", dpsibi, depsbi, dra);
    }
    else {
        return NULL;
    }
}

static PyObject *
_erfa_fk5hip(PyObject *self)
{
    double r5h[3][3], s5h[3];
    eraFk5hip(r5h, s5h);
    return Py_BuildValue("((ddd)(ddd)(ddd))(ddd)",
        r5h[0][0],r5h[0][1],r5h[0][2],
        r5h[1][0],r5h[1][1],r5h[1][2],
        r5h[2][0],r5h[2][1],r5h[2][2],
        s5h[0],s5h[1],s5h[2]);
}

static PyMethodDef _erfaMethods[] = {
    {"bi00", (PyCFunction)_erfa_bi00, METH_NOARGS, "bi00 doc"},
    {"fk5hip", (PyCFunction)_erfa_fk5hip, METH_NOARGS, "fk5hip doc"},
    {NULL, NULL, 0, NULL}
};

static void 
_erfa_d2dtf(char **args, npy_intp *dimensions,
            npy_intp* steps, void* data)
{
    int j, status;
    int in, iy, im, id, ihmsf[4];
    double d1, d2;
    npy_intp i;
    npy_intp N = dimensions[0];
    char *in1 = args[0], *in2 = args[1], *in3 = args[2];
    char *out1 = args[3], *out2 = args[4], *out3 = args[5], 
         *out4 = args[6], *out5 = args[7], *out6 = args[8], *out7 = args[9];
    npy_intp in1_step = steps[0], in2_step = steps[1], in3_step = steps[2];
    npy_intp out1_step = steps[3], out2_step = steps[4],
             out3_step = steps[5], out4_step = steps[6],
             out5_step = steps[7], out6_step = steps[8], 
             out7_step = steps[9];
 
    for (i = 0; i < N; i++) {
        in = *(int *)in1;
        d1 = *(double *)in2;
        d2 = *(double *)in3;

        status = eraD2dtf("UTC", in, d1, d2,
			  &iy, &im, &id, ihmsf);
        *(int *)out1 = iy;
        *(int *)out2 = im;
        *(int *)out3 = id;
	*(int *)out4 = ihmsf[0];
        *(int *)out5 = ihmsf[1];
        *(int *)out6 = ihmsf[2];
        *(int *)out7 = ihmsf[3];

        in1 += in1_step;
        in2 += in2_step;
        in3 += in3_step;
        out1 += out1_step;
        out2 += out2_step;
        out3 += out3_step;
        out4 += out4_step;
	out5 += out5_step;
        out6 += out6_step;
        out7 += out7_step;
    }
}

static PyUFuncGenericFunction d2dtf_func[1] = {&_erfa_d2dtf};
static char d2dtf_types[] = {NPY_INT,  NPY_DOUBLE, NPY_DOUBLE,
			     NPY_INT, NPY_INT, NPY_INT,
			     NPY_INT, NPY_INT, NPY_INT, NPY_INT};
static void *d2dtf_data[1] = {(void *)NULL};

'''

function_template = '''
static void 
_erfa_%s(char **args, npy_intp *dimensions, npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp N = dimensions[0];
    /* declare loop args */
%s
    /* declare steps */
%s
    /* declare function args */
%s
    for (i = 0; i < N; i++) {
    /* assign input args */
%s
    /* call erfa function */
%s
    /* check error */
    /* assign output args */
%s
    /* increment steps */
%s
    }
}
static PyUFuncGenericFunction %s_func[1] = {&_erfa_%s};
static char %s_types[] = {%s};
static void *%s_data[] = {(void *)NULL};
'''

dict_init_template = '''
/* functions dictionnary initialization */
static void
init_dict(PyObject *dict)
{
    PyObject *f;
    int num;
    %s
    return;
}
'''

dict_template = '''
    num = sizeof(%s_func) / sizeof(%s_func[0]);
    f = PyUFunc_FromFuncAndData(%s_func, %s_data, %s_types, num,
                                %d, %d, PyUFunc_None, "%s", "%s", 0);
    PyDict_SetItemString(dict, "%s", f);
    Py_DECREF(f);
'''

init_module_template = '''
#if defined(NPY_PY3K)
static struct PyModuleDef _erfamodule = {
    PyModuleDef_HEAD_INIT,
    "_erfa",
    "Python wrapper for ERFA",
    -1,
    _erfaMethods,
    NULL,
    NULL,
    NULL,
    NULL
};
#endif

#if defined(NPY_PY3K)
#define RETVAL m
PyMODINIT_FUNC PyInit__erfa(void)
#else
#define RETVAL
PyMODINIT_FUNC init_erfa(void)
#endif
{
    PyObject *m, *d;
#if defined(NPY_PY3K)
    m = PyModule_Create(&_erfamodule);
#else
    m = Py_InitModule("_erfa", _erfaMethods);
#endif
    if (!m) {
        return RETVAL;
    }
    import_array();
    //import_umath();
    import_ufunc();
    /* may need some error check here */
    d = PyModule_GetDict(m);
    init_dict(d);

    _erfaError = PyErr_NewException("_erfa.error", NULL, NULL);
    Py_INCREF(_erfaError);
    PyModule_AddObject(m, "error", _erfaError);
    PyModule_AddObject(m, "DPI", PyFloat_FromDouble(ERFA_DPI));
    PyModule_AddObject(m, "D2PI", PyFloat_FromDouble(ERFA_D2PI));
    PyModule_AddObject(m, "DR2D", PyFloat_FromDouble(ERFA_DR2D));
    PyModule_AddObject(m, "DD2R", PyFloat_FromDouble(ERFA_DD2R));
    PyModule_AddObject(m, "DR2AS", PyFloat_FromDouble(ERFA_DR2AS));
    PyModule_AddObject(m, "DAS2R", PyFloat_FromDouble(ERFA_DAS2R));
    PyModule_AddObject(m, "DS2R", PyFloat_FromDouble(ERFA_DS2R));
    PyModule_AddObject(m, "TURNAS", PyFloat_FromDouble(ERFA_TURNAS));
    PyModule_AddObject(m, "DMAS2R", PyFloat_FromDouble(ERFA_DMAS2R));
    PyModule_AddObject(m, "DTY", PyFloat_FromDouble(ERFA_DTY));
    PyModule_AddObject(m, "DAYSEC", PyFloat_FromDouble(ERFA_DAYSEC));
    PyModule_AddObject(m, "DJY", PyFloat_FromDouble(ERFA_DJY));
    PyModule_AddObject(m, "DJC", PyFloat_FromDouble(ERFA_DJC));
    PyModule_AddObject(m, "DJM", PyFloat_FromDouble(ERFA_DJM));
    PyModule_AddObject(m, "DJ00", PyFloat_FromDouble(ERFA_DJ00));
    PyModule_AddObject(m, "DJM0", PyFloat_FromDouble(ERFA_DJM0));
    PyModule_AddObject(m, "DJM00", PyFloat_FromDouble(ERFA_DJM00));
    PyModule_AddObject(m, "DJM77", PyFloat_FromDouble(ERFA_DJM77));
    PyModule_AddObject(m, "TTMTAI", PyFloat_FromDouble(ERFA_TTMTAI));
    PyModule_AddObject(m, "DAU", PyFloat_FromDouble(ERFA_DAU));
    PyModule_AddObject(m, "CMPS", PyFloat_FromDouble(ERFA_CMPS));
    PyModule_AddObject(m, "AULT", PyFloat_FromDouble(ERFA_AULT));
    PyModule_AddObject(m, "DC", PyFloat_FromDouble(ERFA_DC));
    PyModule_AddObject(m, "ELG", PyFloat_FromDouble(ERFA_ELG));
    PyModule_AddObject(m, "ELB", PyFloat_FromDouble(ERFA_ELB));
    PyModule_AddObject(m, "TDB0", PyFloat_FromDouble(ERFA_TDB0));
    PyModule_AddObject(m, "SRS", PyFloat_FromDouble(ERFA_SRS));
    PyModule_AddObject(m, "WGS84", PyLong_FromLong(ERFA_WGS84));
    PyModule_AddObject(m, "GRS80", PyLong_FromLong(ERFA_GRS80));
    PyModule_AddObject(m, "WGS72", PyLong_FromLong(ERFA_WGS72));

    /* ... init for struct type ... */
    
    return RETVAL;
}
'''

NPY_TYPES = {'int':'NPY_INT',
             'double':'NPY_DOUBLE',
             'char':'NPY_STRING'}

no_args = []
double = []
void = []
integer = [] 
char_args = []
astrom_args = []
ldbody_args = []

erfa = ErfaWrapp()
erfa_names = list(erfa.functions.keys())
erfa_names.sort()

_ = open('_erfamodule.c', 'w')
_.write(header_templates)
_dict = '''    f = PyUFunc_FromFuncAndData(d2dtf_func, d2dtf_data, d2dtf_types,
                                     1, 3, 7,
                                     PyUFunc_None, "d2dtf",
                                     "d2dtf docstring",
				     0);

    PyDict_SetItemString(dict, "d2dtf", f);
    Py_DECREF(f);

'''

for f in erfa_names:
    if not erfa.functions[f]['args_in']:
        no_args.append(f)
        continue
    for a in erfa.functions[f]['args']:
        if 'char' in a:
            char_args.append(f)
            continue
        if 'eraASTROM' in a:
            astrom_args.append(f)
            continue
        if 'eraLDBODY' in a:
            ldbody_args.append(f)
            continue
    ## skip these functions for now 
    if f in astrom_args+ldbody_args+no_args+char_args+['eraRx','eraCr', 'eraRz','eraRy']:
        continue
    name = f[3:].lower()
    nb_args = 0 
    nb_args_in = 0
    nb_args_out = 0
    args = ''
    _in = ''
    _out = ''
    era = ''
    args_in = ''
    args_out = ''
    step_in = ''
    step_out = ''
    inc_step = ''
    types = ''
    if erfa.functions[f]['returned_value']:
        args += '%s ret;\n'%erfa.functions[f]['type']
        era += 'ret = '
    era += '%s('%f
    for i in erfa.functions[f]['args']:
        args += '%s;\n'%(i.replace('*', ''))
        shape = 1
        t, n = i.split()[-2:]
        na = n
        if '*' in n:
            na = n.replace('*', '&')
        elif '[' in n:
            n = n.split('[')
            n = [j.strip(']') for j in n]
            na = n[0]
            shape = [int(i) for i in n[1:]]
        era += '%s, '%na
        na = na.lstrip('&')
        if na in erfa.functions[f]['args_in']:
            if shape == 1:
                types += '%s,'%NPY_TYPES[t]
                args_in += 'char *in%d = args[%d];\n'%((nb_args_in+1), nb_args)
                step_in += 'npy_intp in%d_step = steps[%d];\n'%((nb_args_in+1), nb_args)
                inc_step += 'in%d += in%d_step;\n'%((nb_args_in+1), (nb_args_in+1))
                _in += '%s = *(%s *)in%d;\n'%(na, t, nb_args_in+1)
                nb_args_in += 1
                nb_args +=1
            else:
                if len(shape) == 1:
                    for k in range(shape[0]):
                        types += '%s,'%NPY_TYPES[t]
                        args_in += 'char *in%d = args[%d];\n'%((nb_args_in+1), nb_args)
                        step_in += 'npy_intp in%d_step = steps[%d];\n'%((nb_args_in+1), nb_args)
                        inc_step += 'in%d += in%d_step;\n'%((nb_args_in+1), (nb_args_in+1))
                        _in += '%s[%d] = *(%s *)in%d;\n'%(na, k, t, nb_args_in+1)
                        nb_args_in += 1
                        nb_args +=1
                elif len(shape) == 2:
                    for l in range(shape[0]):
                        for k in range(shape[1]):
                            types += '%s,'%NPY_TYPES[t]
                            args_in += 'char *in%d = args[%d];\n'%((nb_args_in+1), nb_args)
                            step_in += 'npy_intp in%d_step = steps[%d];\n'%((nb_args_in+1), nb_args)
                            inc_step += 'in%d += in%d_step;\n'%((nb_args_in+1), (nb_args_in+1))
                            _in += '%s[%d][%d] = *(%s *)in%d;\n'%(na, l, k, t, nb_args_in+1)
                            nb_args_in += 1
                            nb_args +=1

        if na in erfa.functions[f]['args_out']:
            if shape == 1:
                types += '%s,'%NPY_TYPES[t]
                args_out += 'char *out%d = args[%d];\n'%((nb_args_out+1), nb_args)
                step_out += 'npy_intp out%d_step = steps[%d];\n'%((nb_args_out+1), (nb_args))
                inc_step += 'out%d += out%d_step;\n'%((nb_args_out+1), (nb_args_out+1))
                _out += '*(%s *)out%d = %s;\n'%(t, nb_args_out+1, na)
                nb_args_out += 1
                nb_args += 1
            else:
                if len(shape) == 1:
                    for k in range(shape[0]):
                        types += '%s,'%NPY_TYPES[t]
                        args_out += 'char *out%d = args[%d];\n'%((nb_args_out+1), nb_args)
                        step_out += 'npy_intp out%d_step = steps[%d];\n'%((nb_args_out+1), nb_args)
                        inc_step += 'out%d += out%d_step;\n'%((nb_args_out+1), (nb_args_out+1))
                        _out += '*(%s *)out%d = %s[%d];\n'%(t, nb_args_out+1, na, k)
                        nb_args_out += 1
                        nb_args += 1
                elif len(shape) == 2:
                    for l in range(shape[0]):
                        for k in range(shape[1]):
                            types += '%s,'%NPY_TYPES[t]
                            args_out += 'char *out%d = args[%d];\n'%((nb_args_out+1), nb_args)
                            step_out += 'npy_intp out%d_step = steps[%d];\n'%((nb_args_out+1), nb_args)
                            inc_step += 'out%d += out%d_step;\n'%((nb_args_out+1), (nb_args_out+1))
                            _out += '*(%s *)out%d = %s[%d][%d];\n'%(t, nb_args_out+1, na, l, k)
                            nb_args_out += 1
                            nb_args += 1
    era = era[:-2] + ');\n'
    if nb_args_out == 0:
        ## arg_out is the function returned value
        args_out = 'char *out1 = args[%d];\n'%nb_args
        step_out += 'npy_intp out1_step = steps[%d];\n'%nb_args
        inc_step += 'out1 += out1_step;\n'
        _out += '*(%s *)out1 = ret;\n'%erfa.functions[f]['type']
        types += '%s,'%NPY_TYPES[erfa.functions[f]['type']]
        nb_args_out += 1
        nb_args += 1
    types = types[:-1]
    _.write(function_template%(name,
                             args_in+args_out,
                             step_in+step_out,
                             args,
                             _in,
                             era,
                             _out,
                             inc_step,
                             name, name,
                             name, types,
                             name))
    _dict += dict_template%(name, name,
                         name, name, name,
                         nb_args_in, nb_args_out,
                         name, "%s doc"%name,
                         name)

_.write(dict_init_template%_dict)
_.write(init_module_template)
_.close()
