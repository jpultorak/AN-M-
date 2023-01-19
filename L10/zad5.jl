using Polynomials


function IntegratePol(n, k, T = Float64)
    roots = [T(i) for i in 0:n if i != k]
    P = fromroots(roots)
    P_integrated = integrate(P)
    return P_integrated(n) - P_integrated(0)
end

function Ak_n(n, k, a, b, T = Float64)
    h = T((b-a)/n)
    sgn = 1
    if mod(n-k, 2) == 1
        sgn = -1
    end
    fact = T(factorial(k))*T(factorial(n-k))
    return sgn * h * (1/fact) *IntegratePol(n, k, BigFloat)
end

function f(x, T = Float64)
    return T(1)/(T(1)+T(x^2))
end

function Newton_Cotes(n, f, a, b, T = Float64)
    h = (T(b)-T(a))/T(n)
    xs = [T(a) + T(k)*h for k in 0:n]
    res = T(0)
    for k in 0:n
        res += Ak_n(n, k, a, b, BigFloat)*f(xs[k+1], BigFloat)
    end 
    return res
end

for n in 2:2:16
    println(n, ": ", Newton_Cotes(n, f, -4, 4, BigFloat))
end
