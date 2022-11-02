## Exercício 3
> Problema No. 3  - Advecção com difusão e forçante 
> Apresentação dos resultados em **(12/10/2022)**

Vamos simular uma situação mais realista do Problema 2 (advecção com fonte do poluente) ao introduzir o efeito da difusão. Portanto, a equação que governa este problema é dada por:

$\frac{\partial C}{\partial t} + U \frac{\partial C}{\partial x} = K \frac{\partial^2C}{\partial x^2} + F$

Use a mesma fonte periódica do Ex 2, mesma malha no espaço e determine o K (kappa) de forma que o tempo de decaimento seja da ordem de 3 horas. Inicialmente coloque o F no tempo n. Faça o 

$\frac{\partial C}{\partial t} \rarr \frac{(C^{n+1} - C^{n-1})}{2 \Delta t}$

Mantenha a CF radiacional do Ex 2. 

- [x] a. Resolva numericamente com o esquema Leapfrog na advecção e a difusão no tempo n-1, forcante no tempo n-1.
- [x] b. Introduza a forçante pelo método splitting
  - [x] 1. Verifique experimentalmente o critério discutido em Doos et al. Figura 8.6  para a estabilidade numérica do esquema, através de variações do U e K (kappa) com F = 0
  - [x] 2. Com o F do Problema 2, discuta o efeito do splitting. Ou seja, compare a solução com o F calculado no tempo n-1 com a solução em 2 passos. No primeiro passo  (\*) calcule somente o efeito da advecção e difusão e no segundo passo calcule a forçante com a estimativa  no primeiro passo (\*).

Apresente o exercício na forma de um paper, ou seja,  com uma 
- [x] introdução sobre o tema, 
- [x] descrição da metodologia, 
- [x] resultados e 
- [x] discussão dos resultados.  Pode retirar toda a discussão teórica  sobre a advecção e forçante. Foque na difusão e na questão da forma de introduzir a forçante.