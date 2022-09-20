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

     



