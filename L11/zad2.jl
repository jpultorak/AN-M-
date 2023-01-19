
function f(x)
    return sin(x)
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
# zad2 f(x) = sin(x)
total = 15
comp_trapezoidal = [composite_Trapezoidal(f, 2^k, 0.0, pi) for k in 0:total]
comp_simpson = [(4*comp_trapezoidal[k+1]-comp_trapezoidal[k])/3 for k in 1:total]

for k in 0:total
    println("comp_trapezoidal (n = ", 2^k, ") ", comp_trapezoidal[k+1], "   ",  abs(2-comp_trapezoidal[k+1]))
end
for k in 1:total
    println("comp_simpson (n = ", 2^k, ") ", comp_simpson[k], "    ",  abs(2-comp_simpson[k]))
end
# zad 11 f(x) sin(x)e^(-x)
function g(x)
    return cos(x)*cos(x)*exp(-x)
end
# comp_trapezoidal = [composite_Trapezoidal(g, 2^k, 0.0, 100) for k in 0:15]
# comp_simpson = [(4*comp_trapezoidal[k+1]-comp_trapezoidal[k])/3 for k in 1:15]


# println("Trapezoidal Composite: ",comp_trapezoidal)
# println("Simpson Composite: ",comp_simpson)

