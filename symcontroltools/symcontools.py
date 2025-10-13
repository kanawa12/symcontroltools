
import sympy as sp

def testprint(tst):
    print(tst)

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

