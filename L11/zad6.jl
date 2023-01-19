function ClenshawCurtis(f,n)
    # Chebyshev extreme points
    x = cos.(pi*(0:n)/n)
  
    fx = f.(x)/(2n)
    println(fx)
    println(fx[n:-1:2])
    println(vcat(fx,fx[n:-1:2]))
    #vcat fx[n:-1:2] bierze punkty 1:n (czyli bez 0 i n+1)

    # # Fast Fourier transform
    # g = real(FFTW.fft(vcat(fx,fx[n:-1:2])))
    # # Chebyshev coefficients
    # a = vcat( g[1], g[2:n]+g[2*n:-1:n+2], g[n+1] )
    # w = zeros(length(a))
    # w[1:2:end] = 2 ./ (1 .- (0:2:n) .^ 2 )
    # LinearAlgebra.dot(w,a)
end

function f(x)
    return x
end
ClenshawCurtis(f, 5)