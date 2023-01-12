using Plots
using HTTP


function show_f()
    xs = 0.0 : 0.005 : 1
    ys = [f(x) for x in xs]
    scatter(xs, ys)
end

function f(x :: Float64)
    r = HTTP.request("GET", "http://roxy.pythonanywhere.com/f3?x=$x")
    return parse(Float64, String(r.body))
end

# Newton-Cotes for n = 1
function Trapezoidal(f, a, b)
    return (b-a)*(f(a)+f(b))/2
end

# Newton-Cotes for n = 2
function Simpson(f, a, b)
    mid = (a+b)/2
    return (b-a)*(f(a)+ 4*f(mid) + f(b))/6
end

function composite_Trapezoidal(f, n, a, b)
    h = (b-a)/n
    sum = 0
    for i in 1:n-1
        sum += f(a + i*h)
    end
    return h*(f(a) + 2*sum + f(b))/2
end

comp_trapezoidal = [composite_Trapezoidal(f, 2^k, 0.0, 1.0) for k in 0:5]
comp_simpson = [(4*comp_trapezoidal[k+1]-comp_trapezoidal[k])/3 for k in 1:5]

println("Trapezoidal: ", Trapezoidal(f, 0.0, 1.0))
println("Simpson: ",Simpson(f, 0.0, 1.0))


println("Trapezoidal Composite: ",comp_trapezoidal)
println("Simpson Composite: ",comp_simpson)


#error
for k in 1:5
    println(1.61*100/(2.88 * 2^4k))
end
#show_f()