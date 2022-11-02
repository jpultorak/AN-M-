using Printf
using Polynomials

function solve_quadratic_equation(a::BigFloat,b::BigFloat,c::BigFloat)
    Δ=b*b-4.0*a*c;
    x1,x2 = 0.0, 0.0;
    
    if (Δ>0.0)
        sΔ = sqrt(Δ)
        if (b>0.0)
            x1 = (-b-sΔ)/(2.0*a)
            x2 = c/x1
        else
            x2 = (-b+sΔ)/(2.0*a)
            x1 = c/x2
        end
    elseif (Δ<0.0)
        sΔ = sqrt(-Δ)
        if (b>0.0)
            x1 = (-b-sΔ*im)/(2.0*a)
            x2 = c/x1
        else
            x2 = (-b+sΔ*im)/(2.0*a)
            x1 = c/x2
        end
    else
        x1 = -b/(2.0*a);
        x2 = -b/(2.0*a);
    end
    return x1,x2;
end     

function Bairstow(a, iters, u0, v0)
    n = size(a)[1]

    roots = zeros(Complex{BigFloat}, n-1)
    root_cnt = 1
    while n > 2
        b = zeros(BigFloat, n)
        c = zeros(BigFloat, n)
        b[n] = a[n]
        
        c[n] = 0
        c[n-1] = a[n]
        u = u0
        v = v0 
        for _ = 1:iters
            b[n-1] = a[n-1] + u*b[n]
            for k = n-2:-1:1
                b[k] = a[k]+ u*b[k+1] + v*b[k+2]
                c[k] = b[k+1] + u*c[k+1] + v*c[k+2]
            end
        
            detJ = c[1]*c[3]-c[2]*c[2]
            Δu = (c[2]*b[2]-c[3]*b[1])/detJ
            Δv = (c[2]*b[1]-c[1]*b[2])/detJ
            u = u + Δu
            v = v + Δv
        end
        a = b[3:end]
        x1, x2 = solve_quadratic_equation(BigFloat(1), -u, -v)
        roots[root_cnt] = x1;
        roots[root_cnt+1] = x2;
        root_cnt = root_cnt + 2
        n = n-2
    end
    
    if n == 2
        println(a)
        roots[root_cnt] = -a[1]/a[2]
    end
    return roots
end

a = [1., 2., 3., 4., 5.]
#a = [-4, 0, 1]
a = map(BigFloat, a)
w_res = Bairstow(a, 10, 0.1, 0.1)

function w(x)
    1+2*x+3*x*x+4*x*x*x+5*x*x*x*x
end


for x in w_res   
    println(x, "\n")
end

w1 = Polynomial([1., 2., 3., 4., 5.])
println(roots(w1))