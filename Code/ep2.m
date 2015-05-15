
A = [1 1 2 1 0 0 0;
     0 1 6 0 1 0 0;
     1 0 0 0 0 1 0;
     0 1 0 0 0 0 1]
b = [8 12 4 6]
c = [1 2 3 4]
m = 4
n = 7
x = [0 0 0 8 12 4 6]

function simplex (A, b, c, m, n, x)

  # criação do vetor de índices
  bind = []
  B = A
  for i = 1:n
    if x(i) != 0
      bind(end + 1) = i
    else 
      B(:, 1) = []
   endif
  endfor
  
  [L, U] = lu(B)
  
  # criação de um vetor com os custos
  v = []
  for j = 1:n
    if bind != j
      y = transpose(A(:, i)) / L
      u = y / U
      v(end + 1) = c(j) - c*transpose(u)
    endif
  endfor
  
  
  
endfunction

[ind v] = simplex(A, b, c, m, n, x)
