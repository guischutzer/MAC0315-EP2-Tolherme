
A = [1 1 2 1 0 0 0;
     0 1 6 0 1 0 0;
     1 0 0 0 0 1 0;
     0 1 0 0 0 0 1];
b = [8 12 4 6];
c = [10 20 30 40 50 60 70];
m = 4;
n = 7;
x = [0 0 0 8 12 4 6];

# função-invólucro da recursão:
# calcula dados que só precisam ser
# mudados a cada iteração
function [ind, v] = simplex(A, b, c, m, n, x)
  
  # criação do vetor de índices - bind;
  # da matriz de bases - B;
  # da parte básica do vetor de custos - cB;
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
  
  # chamada da função recursiva (iteração == 0)
  [ind v] = simplex_rec(A, b, c, m, n, x, B, bind, cB, 0);
  
endfunction

function [ind, v] = simplex_rec(A, b, c, m, n, x, B, bind, cB, it)

  # decomposição LU de B
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
  
  # determinação dos custos reduzidos
  for j = 1:n
    # caso j não seja variável básica
    if bind != j
      # encontrar p t. q. c_B * B^(-1) = p'
      # com a decomposição LU
      y = cB / U;
      p = y / L;
      # "append" no vetor de custos reduzidos
      # cada cj-barra
      cjbarra(end + 1) = c(j) - p * A(:, j);
      
      # determina a j-ésima variável (não básica)
      # a ser tomada como direção de redução de custos:
      # - é a primeira variável que corresponde a cj-barra
      # negativo;
      # - loop não ocorre caso já esteja estabelecida (j != 0)
      if l == 0 && cjbarra(end) < 0
        l = j;
        # cálculo do vetor u = -dB, com as componentes
        # relacionadas aos índices de mesma posição em bind
        y = transpose(A(:, j)) / transpose(L);
        u = y / transpose(U);
        
        # cálculo do menor valor possível para theta*
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
      # cj-barra correspondente a uma variável básica
      # será 0; passo necessário para organização dos dados
      cjbarra(end + 1) = 0;
        
    endif     
  endfor
  
  # Impressão dos resultados da iteração
  
  printf("\n-------------------------------------\nIterando %d\n\n", it);
  printf("Variáveis Básicas:\n");
  
  for i = sort(bind)
    printf("%d %.5f \n", i, x(i));
  endfor      
    
  printf("\nValor função objetivo: %.5f\n\n", c*transpose(x));
   
  printf("Custos reduzidos:\n")
  for i = 1:n
    if i != bind
      printf("%d %.5f\n", i, cjbarra(i));
    endif
  endfor
  
  if l == 0
    ind = 0;
    v = transpose(x);
    printf("\n=====================================\nSolução ótima encontrada com custo %.5f:\n", c*transpose(x));
    for i = 1:n
      printf("%d %.5f\n", i, x(i))
    endfor
    return
  endif
  
  d = zeros(n, 1);
  for i = 1:m
    d(bind(i)) = -u(i);
    u(i);
  endfor
  d(l) = l;
  
  if startheta == Inf
    ind = -1;
    v = x;
    printf("\n=====================================\nProblema ilimitado. Custo ótimo -∞ na direção:\n");
    for i = 1:n
      printf("%d %.5f\n", i, d(i));
    endfor
    return
  endif
    
  printf("\nEntra na base: %d\n\n", l);
  
  printf("Direção:\n");
  for i = 1:m
    printf("%d %.5f \n", bind(i), -u(i));
  endfor
    
  printf("\nTheta*\n%.5f\n\n", startheta);
    
  printf("Sai da base: %d\n", bind(imin));
    
  # Atualização: novos B, x, cB, bind
  B(:, imin) = A(:, l);
  
  for i = 1:m
    x(bind(i)) = x(bind(i)) - startheta * u(i);
  endfor
  x(l) = startheta;
  
  cB(imin) = c(l);
  
  bind(imin) = l;
    
  it = it + 1;
  
  [ind v] = simplex_rec(A, b, c, m, n, x, B, bind, cB, it);
  

endfunction

[ind v] = simplex(A, b, c, m, n, x);