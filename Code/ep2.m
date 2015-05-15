
A = [1 1 2 1 0 0 0;
       0 1 6 0 1 0 0;
       1 0 0 0 0 1 0;
       0 1 0 0 0 0 1]
b = [8 12 4 6]
c = [10 20 30 40 50 60 70]
m = 4
n = 7
x = [0 0 0 8 12 4 6]

simplex_rec = max_recursion_depth(5);

function [ind, v] = simplex(A, b, c, m, n, x)
  
  # criação do vetor de índices
  bind = [];
  B = A;
  for i = fliplr(1:n);
    if x(i) != 0
      bind(end + 1) = i;
      cB(end + 1) = c(i);
    else 
      B(:, i) = [];
   endif
  endfor
  bind = fliplr(bind);
  cB   = fliplr(cB);
  
  [ind v] = simplex_rec(A, b, c, m, n, x, B, bind, cB, 0);
  
endfunction

function [ind, v] = simplex_rec(A, b, c, m, n, x, B, bind, cB, it)

  [L, U] = lu(B);

  # criação de um vetor com os custos cj-barra
  # associados às variáveis não-básicas
  cjbarra = [];
  # inicialização do theta* em infinito
  startheta = Inf;
  # inicialização da variável a ser escolhida
  # como direção não-básica
  l = 0;
  # inicialização do vetor u = -dB
  u = [];
  imin = 0;
  for j = 1:n
    if bind != j
      # encontrar p' t. q. c'B * B^(-1) = p'
      y = cB / transpose(U);
      p = y / transpose(L);
      cjbarra(end + 1) = c(j) - p * A(:, j);
      
      if l == 0 && cjbarra(end) < 0
        l = j;
        y = transpose(A(:, j)) / L;
        u = y / U;
        
        for i = 1:m
          if u(i) > 0
            if x(bind(i))/u(i) < startheta
              startheta = x(bind(i))/u(i);
              imin = i;
            endif
          endif
        endfor
      endif
      
      else 
        cjbarra(end + 1) = 0;
        
    endif
      
  endfor
    
  printf("----------------\nIterando %d\n\n", it);
  printf("Variaveis Basicas:\n");
  for i = bind
    printf("%d %.5f \n", i, x(i));
  endfor      
    
  printf("\nValor funcao objetivo: %.5f\n\n", c*transpose(x));
   
  printf("Custos reduzidos:\n")
  for i = 1:n
    if i != bind
      printf("%d %.5f\n", i, cjbarra(i));
    endif
  endfor
  
  if l == 0
    ind = 0;
    v = x;
    return
  endif
  
  if startheta == Inf
    ind = -1;
    v = x;
    return
  endif
    
  printf("\nEntra na base: %d\n\n", l);
    
  printf("Direcao:\n");
  for i = 1:m
    printf("%d %.5f \n", bind(i), -u(i));
  endfor
    
  printf("\nTheta*\n%.5f\n\n", startheta);
    
  printf("Sai da base: %d\n", bind(imin));
    
  
  B(:, imin) = A(:, l)
  
  for i = 1:m
    x(bind(i)) = x(bind(i)) - startheta * u(i);
  endfor
  x(l) = startheta;
  
  cB(imin) = c(l);
  
  bind(imin) = l;
    
  it = it + 1;
  
  [ind v] = simplex_rec(A, b, c, m, n, x, B, bind, cB, it);
  return
  

endfunction

[ind v] = simplex(A, b, c, m, n, x);
