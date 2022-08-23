# Modelagem Atmosférica - 2o semestre 2022
Exercicíos desenvolvidos em Python como parte dos estudos de doutorado em meteorologia.

## Exercicío 1
- [ ] Resolver problema de transporte de um pulso inicial de poluição numa situação em que o vento seja constante.
   - [ ] Desenvolver as soluções numéricas: Malha de espaçamento, condição inicial, U = 10 m/s. Integrar durante o tempo suficiente para que perturbação inicial volte para a parte central da malha no caso da fronteira; demais casos, integrar até que a perturbação chegue totalmente na fronteira.
   - [ ] Método de solução
     1. Inicie com um esquema de ordem 1 (cuidado com a estabilidade do esquema na escolha do delta t na derivada espacial).
     2. Compare a solução dada em 1 com um esquema leapfrog de ordem 2 (espaço e tempo)
     3. Modifique a derivada no espaço para um esquema de 4a ordem (será necessário manter a segunda ordem no ponto imediatamente vizinho da fronteira no caso da CF fixa).
     4. Use também um método iterativo, por exemplo, o `Esquema Matsun`.
     5. Use um `método implícito`.
     6. Finalmente, use um Runge Kutta de 4ª ordem no tempo com um esquema de 4[ ordem no espaço
     
     Obs: quanto a CF, não é necessário fazer os 3 casos para todos os experimentos. Faça apenas para o caso Leapfrog. Nos demais, faça sempre a CF periódica.
     
- [ ] Apresentar o exercício na forma de um "paper": introdução, metodologia, resultados e discussão dos resultados.
     



