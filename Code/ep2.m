
  A = [1 1 2 1 0 0 0;
       0 1 6 0 1 0 0;
       1 0 0 0 0 1 0;
       0 1 0 0 0 0 1]
  b = [8 12 4 6]
  c = [10 20 30 40 50 60 70]
  m = 4
  n = 7
  x = [0 0 0 8 12 4 6]

  function retval = simplex (A, b, c, m, n, x, it)
  
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
    # inicialização do theta* em infinito
    startheta = Inf
    # inicialização da variável a ser escolhida
    # como direção não-básica
    l = 0
    # inicialização do vetor u = -dB
    u = []
    imin = 0
    for j = 1:n
      if bind != j
        # encontrar p' t. q. c'B * B^(-1) = p'
        y = cB / transpose(U)
        p = y / transpose(L)
        cjbarra(end + 1) = c(j) - p * A(:, j)
        
        if l == 0 && cjbarra(end) < 0
          l = j
          y = transpose(A(:, j)) / L
          u = y / U
          
          for i = 1:m
            if u(i) > 0
              if x(bind(i))/u(i) < startheta
                startheta = x(bind(i))/u(i)
                imin = bind(i)
              endif
            endif
          endfor
        endif
        
        else 
          cjbarra(end + 1) = 0
        
      endif
      
    endfor
    
    #if l == 0
    #  retval = 0
    #  return
    #endif
    
    #if startheta == Inf
    #  retval = -1
    #  return
    #endif
    
    printf("Iterando %d\n\n", it)
    printf("Variaveis Basicas:\n")
    for i = bind
      printf("%d %.5f \n", i, x(i))
    endfor      
    
    printf("\nValor funcao objetivo: %.5f\n\n", c*transpose(x))
    
    printf("Custos reduzidos:\n")
    for i= 1:n
      if i != bind
        printf("%d %.5f\n", i, cjbarra(i))
      endif
    endfor
    
    printf("\nEntra na base: %d\n\n", l)
    
    printf("Direcao:\n")
    for i = 1:m
      printf("%d %.5f \n", i, -u(i))
    endfor
    
    printf("\nTheta*\n %.5f\n\n", startheta)
    
    printf("Sai da base: %d\n", imin)
    
    
    
      
    
    retval = 1
    
  endfunction

  [v] = simplex(A, b, c, m, n, x, 0)
