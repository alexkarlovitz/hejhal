import sympy as sp

if __name__ == '__main__' :
    x, y = sp.symbols('x y')
    r = sp.sqrt( (x**2 + (y-1)**2) / (x**2 + (y+1)**2) )

    r_x = sp.diff(r, x)
    r_xx = sp.diff(r_x, x)
    r_y = sp.diff(r, y)
    r_yy = sp.diff(r_y, y)

    sp.init_printing()
    coeff = r_xx + r_yy
    print sp.simplify(coeff)
