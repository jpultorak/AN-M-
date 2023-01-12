using Plots, Polynomials

function f(x)
    @assert -1 <= x <=1 "point outside domain"
    return 1/(25*x^2 + 1)
end

function max_error(w, N=1000)
    mx = 0
    for i in 0:1000
        x_k = -1 + 2*i/N
        #println(x_k, "      ", f(x_k), "     ", w(x_k))
        val = abs(f(x_k)-w(x_k))
        mx = max(mx, val)
    end
    return mx
end

function interpolate(xs)
    ys = map(f, xs)
    #println(xs)
    #println(ys)
    P = fit(xs, ys)
    #println(P)
    return max_error(P)
end


xs_0 = [-1 + 2*x/9 for x in 0:9]
T9_extrema = [cos(i*pi/9) for i in 0:9]
T10_zeros = [cos(i*pi/10) for i in 1:10]
#println(T9_extrema)


Ta = Polynomial([1])
Tb = Polynomial([0, 1])
for i in 2:9
    Tc = Polynomial([0, 2])*Tb - Ta
    global Ta= Tb
    global Tb= Tc
    println(Tc)
end
#println(interpolate(xs_0))
#println(interpolate(T9_extrema))
#println(interpolate(T10_zeros))