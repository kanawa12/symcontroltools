
import sympy as sp

def testprint(tst):
    print(tst)


def dp(ep, capt="", pmode=0, isptype = False):
	"""
	depend :
	import sympy as sp
	import os
	"""
	pstr = capt +": "
	if isptype:
		pstr += str(type(ep))
		if hasattr(ep, "shape"):
			pstr += ", " + str(ep.shape)
	dstr = "$$" + sp.latex(ep) + "$$" + os.linesep + "<details><summary>sympy.srepr()</summary><p>" + sp.srepr(ep) + "</p></details>"
	if capt != "": 
		print(pstr)
	display(ep); print(dstr)
	return


def setf(tp, suf, size):
    """
    連番の添え字付きのシンボリック変数を生成するテキストを生成する。
    tpは本体の文字、
    例: {a_{p}}_{0}, ..., {a_{p}}_{2} : setf("a", "p", 3)
    """
    return "{" + tp + "_{"+suf + "}}" + "_{1:"+str(size + 1)+"}"

def get_SISO_sims(dim, suf):
    A = sp.Matrix( a := sp.symbols(setf("a", suf, dim**2)) ).reshape(dim, dim)
    B = sp.Matrix( b := sp.symbols(setf("b", suf, dim)) ).reshape(dim, 1)
    C = sp.Matrix( c := sp.symbols(setf("c", suf, dim)) ).reshape(1, dim)
    D = sp.Matrix( [0] )
    return [(a, b, c), (A, B, C, D)]


def eqs_to_mat(eqList: list, valList: list) -> sp.Matrix:
	"""
	シンボリック式のリストを行列に変換する。
	ex:
	eqList : List
	valList = [ap[1], ap[3], bp[0], bp[1]]
	"""
	valList.append(1)
	return sp.Matrix([sp.Matrix([ sp.Poly(eq, valList[:-1]).coeff_monomial(v) for v in valList ]).reshape(1,len(valList)) for eq in list(eqList) ])

def get_tf(A, B, C, s):
    tf = C * (s*sp.eye(A.shape[0]) - A).inv() * B
    return tf
