
A = [1 1 2 1 0 0 0;
     0 1 6 0 1 0 0;
     1 0 0 0 0 1 0;
     0 1 0 0 0 0 1]
b = [8 12 4 6]
c = [10 20 30 40 50 60 70]
m = 4
n = 7
x = [0 0 0 8 12 4 6]

function retval = simplex (A, b, c, m, n, x)

  # criação do vetor de índices
  bind = []
  B = A
  for i = fliplr(1:n)
    if x(i) != 0
      bind(end + 1) = i
      cB(end + 1) = c(i)
    else 
      B(:, i) = []
   endif
  endfor
  bind = fliplr(bind)
  cB   = fliplr(cB)
  
  [L, U] = lu(B)
  
  # criação de um vetor com os custos cj-barra
  # associados às variáveis não-básicas
  cjbarra = []
  for j = 1:n
    if bind != j
      # encontrar p' t. q. c'B * B^(-1) = p'
      y = cB / transpose(U)
      p = y / transpose(L)
      cjbarra(end + 1) = c(j) - p * A(:, j)
    else 
      cjbarra(end + 1) = 0
    endif
  endfor
  
  if cjbarra >= 0
    return
  endif
  
  
  retval = 1
  
endfunction

[v] = simplex(A, b, c, m, n, x)
