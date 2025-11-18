import sympy as sp
import os
import numpy as np
import pandas as pd
import control as ct
import io
import json
from pathlib import Path


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
	#dstr = "$$" + sp.latex(ep) + "$$" + os.linesep + "<details><summary>sympy.srepr()</summary><p>" + sp.srepr(ep) + "</p></details>"
	dstr = "$$" + sp.latex(ep) + "$$" + "<details><summary>sympy.srepr()</summary><p>" + sp.srepr(ep) + "</p></details>"
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


def eqs_to_mat(eqList: list, valList: list) -> sp.Matrix:
	"""
	シンボリック式のリストを行列に変換する。
	ex:
	eqList : List
	valList = [ap[1], ap[3], bp[0], bp[1], 1]
	"""
	return sp.Matrix([sp.Matrix([ sp.Poly(eq, valList[:-1]).coeff_monomial(v) for v in valList ]).reshape(1,len(valList)) for eq in list(eqList) ])

def eqs_to_mateqs(eqList: list, valList: list) -> sp.Matrix:
	"""
	aaa
	"""
	mat = eqs_to_mat(eqList, valList)
	matA = mat[:, :-1]
	matX = sp.Matrix(valList[:-1])
	matY = -1*mat[:, -1]
	return ( sp.Eq(sp.MatMul(matA, matX, evaluate=False), matY, evaluate=False), matA, matX, matY )


def getDataDict(csv_order, csv_data: str, isFile: bool):
	if isFile:
		df = pd.read_csv(csv_data, header=None, names=csv_order, skipinitialspace=True, dtype=float)
	else:
		f = io.StringIO(csv_data.strip())
		df = pd.read_csv(f, header=None, names=csv_order, skipinitialspace=True, dtype = float)
	datadict = df.to_dict(orient='list')
	return datadict

def setDataOrder(arg_order, datadict):
	ordered_values = [datadict[key] for key in arg_order]
	val = [ list(d) for d in list(zip(*ordered_values))]
	return val

def get_tf(A, B, C, s):
    tf = C * (s*sp.eye(A.shape[0]) - A).inv() * B
    return tf


def ss_to_ssobs(sys):
	n=sys.A.shape[0]
	T = np.zeros([n,n])
	for i in range(n):
	    T[i, (n-1)-i] = 1
	return ct.similarity_transform(sys, T)

def zpk_to_ss_obs(zeros, poles, gein):
	zpksys = ct.zpk(zeros=zeros, poles=poles, gain=gain)
	sssys  = ct.tf2ss(zpksys)

def sym_to_ssobs(tfs, sym):
	nc = np.array(sp.poly(sp.numer(tfs), sym).all_coeffs())
	dc = np.array(sp.poly(sp.denom(tfs), sym).all_coeffs())
	nc = [float(i) for i in nc]
	dc = [float(i) for i in dc]
	edc = ct.tf(*[nc, dc])
	sssys  = ct.tf2ss(edc)
	return ss_to_ssobs(sssys)

def readsacjson(keys):
	script_path = Path(__file__).resolve()
	script_dir = script_path.parent
	json_path = script_dir / "srepr_texts.json"
	with open(json_path, 'r', encoding='utf-8') as f:
		d = json.load(f)
	return [d[k] for k in keys]
	











