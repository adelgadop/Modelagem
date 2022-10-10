# Modelagem Atmosférica - 2o semestre 2022
Exercicíos desenvolvidos em Python como parte dos estudos de doutorado em meteorologia.
## Exercicío 1
- [x] Resolver problema de transporte de um pulso inicial de poluição numa situação em que o vento seja constante.
   - [x] Desenvolver as soluções numéricas: Malha de espaçamento, condição inicial, U = 10 m/s. Integrar durante o tempo suficiente para que perturbação inicial volte para a parte central da malha no caso da fronteira periódica; demais casos, integrar até que a perturbação chegue totalmente na fronteira.
   - [x] Solução analítica para função Gaussiana
   - [x] Solução analítica para função retângulo.
   - [x] Método de solução
     1. [x] Inicie com um esquema de ordem 1 (cuidado com a estabilidade do esquema na escolha do delta t na derivada espacial).
        1. [x] Periódica (ok)
     2. [x] Compare a solução dada em 1 com um esquema leapfrog de ordem 2 (espaço e tempo)
        1. [x] Fixa (ok)
        2. [x] Periódica (ok)
        3. [x] Radiacional (ok)
     3. [x] Modifique a derivada no espaço para um esquema de 4a ordem (será necessário manter a segunda ordem no ponto imediatamente vizinho da fronteira no caso da CF fixa).
        1. [x] Periódica
     4. [x] Use também um método iterativo, por exemplo, o `Esquema Matsuno`.
        1. [x] Periódica
     5. [x] Use um `método implícito` (e.g., Crank-Nicolson scheme).
        1. [x] Periódica
     6. [x] Finalmente, use um Runge Kutta de 4ª ordem no tempo com um esquema de 4to ordem no espaço
        1. [x] Periódica
     
     Obs: quanto a CF, não é necessário fazer os 3 casos (fixa, radiacional e periódica) para todos os experimentos. Faça apenas para o caso Leapfrog. Nos demais, faça sempre a CF periódica.
     
- [x] Apresentar o exercício na forma de um "paper": introdução, metodologia, resultados e discussão dos resultados.    
  - Para gerar o pdf usando pandoc:
   `pandoc -V lang=pt -H format.sty paper.md --filter pandoc-crossref --citeproc -o paper.pdf`
  - [x] Introdução
  - [x] Metodologia
  - [x] Resultados e discussão
  - [x] Conclussões preliminares
---
## Exercicío 2
Vamos supor que você tenha uma fonte de poluição pontual (pe. uma chaminé) que emite em forma de pulso senoidal (definido abaixo) com o mesmo campo básico U = 10 $m.s^{-1}$ do exercício 1 (eq. da advecção linearizada). Por tanto, a equação que governa este problema é dado por:

$\frac{\partial C}{\partial t}+U \frac{\partial C}{\partial x}=F$ onde $F=F(x,t)$

$F(i=100, n\Delta t)= sin(\omega .n\Delta t)$ para n = 0, ... N$_{max}$ se $sin(\omega .n.\Delta t) > 0$

$F(i=100, n\Delta t)=0$ caso $sin(\omega .n\Delta t)< 0$ como $\omega = \frac{2\pi}{1800}$ $s^{-1}$
- [x] Discuta a forma da solução analítica no diagrama de Hovmoller (x,t)
- [x] Resolva numéricamente com os esquemas:
  - [x] (a) Avançado no tempo atrasado no espaço
  - [x] (b) Leapfrog 2ª ordem
  - [x] (c) 4ª ordem no espaço
  - [x] (d) Implícito (e.g., Crank-Nicolson)
  
com a mesma malha do `Exercício 1`. Porém, use uma CI identicamente nula em todo o domínio e com a fonte da grandeza $F$ dada pela expressão acima no ponto central da malha.
Como condição de fronteira do lado direito da malha, no caso Leapfrog e 4ª ordem, a solução dada pelo esquema de 1ª ordem (a). Observe que este esquema permite que a perturbação passe pela fronteira, sem reflexão.
- [x] (e) É realista a solução no caso Leapfrog e 4ª ordem? Explique por qué a concentração não é nula à esquerda da fonte. Como tornar o resultado mais realista? (dica... filtragem...)
- [x] (f) O que acontece com o método implícito? Observe que é incondicionalmente estável. Portanto, o $\Delta t$ pode tornar o CFL > 1 e manter a estabilidade. O que acontece com a solução a medida que $\Delta t$ aumenta acima de 1 (experimente com CFL=2, 4...)

- [x] Apresente o exercício na forma de um "paper", ou seja:
  - [x] Introdução sobre o tema
  - [x] Descrição da metodologia,
  - [x] Resultados
  - [x] Discussão dos resultados. 

Para gerar o pdf usando pandoc:
   `pandoc -H format.sty paper.md --filter pandoc-crossref --citeproc -o paper.pdf`
#### Dica
Explorar modo computacional no espaço, liga a fonte o vento vai transportar na direção. Vai aparecer concentrações negativas na retaguarda. A ideia é explorar os filtros no espaço. Alguns esquemas precisam filtros no espaço e no tempo.

**Pegadinha**: A fonte tem uma oscilação fonte é positiva. Ficar atento no valor 1800, que impacta o CFL, não vai representar a fonte.

**Brincadeira**: Cada 1800 vai emitir a chaminé. O objetivo é $\Delta t$ pensar no fenómeno. Cuidado com o esquema implícito na relação $\Delta t$ e CFL.

## Exercicío 3
> Problema No. 3  - Advecção com difusão e forçante
> 
> Apresentação dos resultados em **(12/10/2022)**

Vamos simular uma situação mais realista do Problema 2 (advecção com fonte do poluente) ao introduzir o efeito da difusão. Portanto, a equação que governa este problema é dada por:

$\frac{\partial C}{\partial t} + U \frac{\partial C}{\partial x} = K \frac{\partial^2C}{\partial x^2} + F$

Use a mesma fonte periódica do Ex 2, mesma malha no espaço e determine o K de forma que o tempo de decaimento seja da ordem de 3 horas. Inicialmente coloque o F no tempo n. Faça o 

$\frac{\partial C}{\partial t} \rarr \frac{(C^{n+1} - C^{n-1})}{2 \Delta t} $

Mantenha a CF radiacional do Ex 2. 

- [ ] a. Resolva numericamente com o esquema Leapfrog na advecção e a difusão no tempo n-1, forcante no tempo n-1.
- [ ] b. Introduza a forçante pelo método splitting
  - [ ] 1. Verifique experimentalmente o critério discutido em Doos et al. Figura 8.6  para a estabilidade numérica do esquema, através de variações do U e K (kapa) com F = 0
  - [ ] 2.	Com o F do Problema 2, discuta o efeito do splitting. Ou seja, compare a solução com o F calculado no tempo n-1 com a solução em 2 passos. No primeiro passo  (\*) calcule somente o efeito da advecção e difusão e no segundo passo calcule a forçante com a estimativa  no primeiro passo (\*).

Apresente o exercício na forma de um paper, ou seja,  com uma 
- [ ] introdução sobre o tema, 
- [ ] descrição da metodologia, 
- [ ] resultados e 
- [ ] discussão dos resultados.  Pode retirar toda a discussão teórica  sobre a advecção e forçante. Foque na difusão e na questão da forma de introduzir a forçante.



