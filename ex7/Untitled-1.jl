using SymPy
function x_1()
    @vars L
    @vars x
    @vars E1
    @vars E2
    @vars E3
    @vars t
    @vars ħ
    func = exp(-im*E1*t/ħ)*sin(π*x/L) + exp(-im*E2*t/ħ)*cos(2π*x/L) + exp(-im*E3*t/ħ)*sin(3π*x/L)
    get = integrate(func*x*func,(x,0,L))
    display(simplify(N(get)))
end # end function


x_1()



