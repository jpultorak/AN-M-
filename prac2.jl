function gaussian_elimination(A)

    if size(A)[1] != size(A)[2]
        println("Matrix must be square")
        return -1
    end

    n = size(A)[1]

    for i in 1:n-1
        if A[i, i] == 0
            println("Gaussian elimination failed")
            return -1
        end

        for j in i+1:n
            l_ij = A[j, i]/A[i, i]
            A[j, i] = l_ij

            for k in i+1:n
                A[j, k] -= l_ij*A[i, k]
            end
        end
    end
    perm = [i for i in 1:n]
    return A, perm
end     

# eliminacja gaussa z skalowanym wyborem wierszy głownych
function gaussian_elimination_with_partial_pivoting(A)

    if size(A)[1] != size(A)[2]
        println("Matrix must be square")
        return -1
    end

    n = size(A)[1]

    perm = [i for i in 1:n]
    s = [maximum(map(abs, A)[i, :]) for i in 1:n]

    for i in 1:n-1
  
        # wybór wiersza głównego (pivoting)
        cur_row = i
        for k in i:n
             if abs(A[perm[k], i]/s[perm[k]]) > abs(A[cur_row, i]/s[cur_row])
                cur_row = k
             end
        end
        # zamiana wierszy [cur_row] oraz [i]
        temp = perm[cur_row]
        perm[cur_row] = perm[i]
        perm[i] = temp
    
        # odejmujemy wiersz p[i] od następnych wierszy
        # pamiętamy mnożniki w zerowanych elementach macierzy
        for j in i+1:n
            l_ij = A[perm[j], i]/A[perm[i], i]
            # zamiast zerować A[perm[j], i] zapisujemy w nim mnożnik l_ij]
            A[perm[j], i] = l_ij
            for k in i+1:n
                A[perm[j], k] -= l_ij*A[perm[i], k]   
            end
        end

    end

    return A, perm
end   

# funkcja rozwiazująca system Ax = b mając rokzład PLU
# podstawianie w tył oraz w przód
function solve_system(A, perm, b)

    n = size(A)[1]    
    # rozwiązanie układu Ly = Pb

    for k in 1:n-1
        for i in k+1:n
            b[perm[i]] -= A[perm[i], k]*b[perm[k]]
        end
    end

    x = zeros(Float64, n)
    for i in n:-1:1
        prev = 0
        for j in i+1:n
            prev += A[perm[i], j]*x[j]
        end
        x[i] = (b[perm[i]]-prev)/A[perm[i], i]    
    end
    return x    
end

function inverse_pivoting(A)
    (A, perm) = gaussian_elimination_with_partial_pivoting(A)
    n = size(A)[1]

    # tworzenie macierzy jednostkowej
    I = zeros(Float64, n, n)
    for i in 1:n
        I[i, i] = 1
    end
    
    A_inv = zeros(Float64, n, n)
    for i in 1:n
        A_inv[:, i] = solve_system(A, perm, I[:, i])
    end
    return A_inv
end

function inverse(A, pivoting = true)
    if pivoting
        (A, perm) = gaussian_elimination_with_partial_pivoting(A)
    else 
        (A, perm) = gaussian_elimination(A)
    end
    n = size(A)[1]

    # tworzenie macierzy jednostkowej
    I = zeros(Float64, n, n)
    for i in 1:n
        I[i, i] = 1
    end
    
    A_inv = zeros(Float64, n, n)
    for i in 1:n
        A_inv[:, i] = solve_system(A, perm, I[:, i])
    end
    return A_inv
end

function verify_inv(A, pivoting = true)
    n = size(A)[1]
    I = zeros(Float64, n, n)
    for i in 1:n
        I[i, i] = 1
    end

    A_cpy = copy(A)
    A_inv = inverse(A_cpy, pivoting)

    display(A)
    display(A*A_inv - I)
    display(A*inv(A)- I)
end

function Hilbert_matrix(n)
    H = zeros(Float64, n, n)
    for i in 1:n
        for j in 1:n
            H[i, j] = 1/(i+j-1)
        end
    end
    return H
end

function Pei_matrix(n, d)
    Pei = ones(Float64, n, n)
    for i in 1:n
        Pei[i, i] = d
    end
    return Pei
end

function test1()
    A_test = [[0.0000000000000000001, 1] [1, 1]]

    display(inverse(A_test))
    display(inverse(A_test, false))
    # verify_inv(A_test, false)
    # verify_inv(A_test, true)
end

display(inverse(Pei_matrix(5, 1.1)))
#test1()
# verify_inv(Pei_matrix(10, 0.999999999))
# println("\n\n NO PIVOTING \n \n")
# verify_inv(Pei_matrix(10, 0.999999999), false)
