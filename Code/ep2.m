
A = [1 1 2 1 0 0 0;
     0 1 6 0 1 0 0;
     1 0 0 0 0 1 0;
     0 1 0 0 0 0 1]
b = [8 12 4 6]
c = [1 2 3 4]
m = 4
n = 7
x = [0 0 0 8 1 2 4 6]

function simplex (A, b, c, m, n, x)

  # criação do vetor de índices
  xB = []
  B = A
  for i = 1:n
    if x(i) != 0
      xB(end + 1) = true
    else 
      xB(end + 1) = false
      B(:, 1) = []
   endif
   
   
   
  endfor
  
  
endfunction

[ind v] = simplex(A, b, c, m, n, x)
