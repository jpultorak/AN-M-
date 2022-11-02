using Printf


function w(x, T = Float64)
    x = T(x)
    return x^3-6x^2+3x-T(0.149)
end


# w(x) function rewritten using Horner schema
function w1(x, T = Float64)
    x = T(x)
    return ((x-6)*x + 3)*x-T(0.149)
end

result = BigFloat(-14.636489)

res1 = [w(4.71, Float16), w(4.71, Float32), w(4.71, Float64)]
res2 = [w1(4.71, Float16), w1(4.71, Float32), w1(4.71, Float64)]

println(res1)
println(res2)

println("\nWersja 1")
for x in map(BigFloat, res1)
    println(abs(x- result)/result)
end
println("\nWersja 2")
for x in map(BigFloat, res2)
    println(abs(x- result)/result)
end