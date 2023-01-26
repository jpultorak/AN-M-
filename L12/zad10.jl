using BenchmarkTools

function gaussian_elimination(A)

    if size(A)[1] != size(A)[2]
        println("Matrix must be square")
        return -1
    end

    n = size(A)[1]
    U = A

    L = zeros(Float64, n, n)
    for i in 1:n
        L[i, i] = 1.0
    end

    for i in 1:n-1
        if U[i, i] == 0
            println("Gaussian elimination failed")
            return -1
        end

        for j in i+1:n
            l_ij = U[j, i]/U[i, i]
            U[j, :] -= l_ij*U[i, :]
            L[j, i] = l_ij
        end
    end
    return L, U
end     

# eliminacja  gaussa z częściowym wyborem wierszy głownych
function gaussian_elimination_with_partial_pivoting(A)

    if size(A)[1] != size(A)[2]
        println("Matrix must be square")
        return -1
    end

    n = size(A)[1]

    U = copy(A)
    L = zeros(Float64, n, n)

    perm = [i for i in 1:n]
    s = [maximum(map(abs, A)[i, :]) for i in 1:n]

    for i in 1:n-1
  
        cur_row = i
        for k in i:n
             if abs(U[perm[k], i]/s[perm[k]]) > abs(U[cur_row, i]/s[cur_row])
                cur_row = k
             end
        end

        temp = perm[cur_row]
        perm[cur_row] = perm[i]
        perm[i] = temp
    

        for j in i+1:n
            l_ij = U[perm[j], i]/U[perm[i], i]
            L[perm[j], i] = l_ij
            for k in i:n
                U[perm[j], k] -= l_ij*U[perm[i], k]   
            end
        end

    end

    for i in 1:n
        L[perm[i], i] = 1
    end

    P = zeros(Int8, n, n)
    for i in 1:n
        P[i, perm[i]] = 1
    end

    return P*L, P*U, P
end   
function solve_system1(L, U, b)

    n = size(L)[1]    
    # rozwiązanie układu Ly = Pb
    y = zeros(Float64, n)

    for k in 1:n
        prev = 0
        for i in 1:k-1
            prev += L[k, i]*y[i]
        end
        y[k] = (b[k]-prev)/L[k, k]
    end

    # rozwiązanie układu Ux = y
    x = zeros(Float64, n)
    for k in n:-1:1
        prev = 0
        for i in k+1:n
            prev += U[k, i]*x[i]
        end
        x[k] = (y[k]-prev)/U[k, k]    
    end
    return x    
end
# funkcja rozwiazująca system Ax = b mając rokzład PLU
# podstawianie w tył oraz w przód
function solve_system(L, U, P, b)

    n = size(L)[1]    
    b = P*b

    # rozwiązanie układu Ly = Pb
    y = zeros(Float64, n)

    for k in 1:n
        prev = 0
        for i in 1:k-1
            prev += L[k, i]*y[i]
        end
        y[k] = (b[k]-prev)/L[k, k]
    end

    # rozwiązanie układu Ux = y
    x = zeros(Float64, n)
    for k in n:-1:1
        prev = 0
        for i in k+1:n
            prev += U[k, i]*x[i]
        end
        x[k] = (y[k]-prev)/U[k, k]    
    end
    return x    
end



#A = [[2., 1., 3.] [-2., 3., -1.] [1., -2., -1.]]
n = 100


A = rand(n, n)
b = rand(Float64, n)

#L, U, P = gaussian_elimination_with_partial_pivoting(A)

function solve_LU(A, b)
    L, U = gaussian_elimination(A)
    return solve_system1(L, U,  b)
end

function solve_builtin(A, b)
    return \(A, b)
end

#err = abs.(x1-x2)
# println(maximum(err))