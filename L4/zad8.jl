using Printf
function steffensen(f, x0, eps=1.0e-20, MaxIt = 20)
    itcnt = 0
    xn = x0
    while(abs(f(xn)) > eps && itcnt < MaxIt)
        @printf("%6.1e\n", f(xn));
        f_xn = f(xn)
        x_nxt = xn - (f_xn*f_xn)/(f(xn + f_xn) - f_xn)
        xn = x_nxt
    end
    @printf("%6.1e\n", f(xn));
    return xn
end

e = 2.7182818284
function f(x)
    e^(-x) - sin(x)
end

println(steffensen(f, BigFloat(3.5)), "\n\n")
println(steffensen(f, BigFloat(-3)))