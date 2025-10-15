import sympy as sp
import os


def dp(ep, capt="", pmode=0, isptype = False):
	"""
	sympyで計算した数式をjupyter notebook上で表示する際、latex表現とsrepr()テキストを同時に表示する。
	引数 pmode : 2のときはそのままテキスト出力 
	isptype : 渡された変数の型と(存在すれば)形を出力するか
	depend :
	import sympy as sp
	import os
	"""
	if pmode == 2:
		display(ep)
		return
	pstr = str(capt) +": "
	if isptype:
		pstr += str(type(ep))
		if hasattr(ep, "shape"):
			pstr += ", " + str(ep.shape)
	dstr = "$$" + sp.latex(ep) + "$$" + os.linesep + "<details><summary>sympy.srepr()</summary><p>" + sp.srepr(ep) + "</p></details>"
	if str(capt) != "": 
		display(pstr)
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



def get_tf(A, B, C, s):
    tf = C * (s*sp.eye(A.shape[0]) - A).inv() * B
    return tf





