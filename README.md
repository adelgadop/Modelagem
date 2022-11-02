# Modelagem Atmosférica - 2o semestre 2022
Exercícios desenvolvidos em Python como parte dos estudos de doutorado em meteorologia:

- [x] **Exercício 1**: Problema de transporte de um pulso inicial numa situação de vento constante com solução analítica (CFL=1) e numérica com diferentes esquemas de aproximação (e.g., ordem 1, leapfrog, Matsuno, Crank-Nicolson e Runge Kutta 4ª ordem). Os experimentos consideraram condições de fronteira tipo fixa, radiacional e periódica.
- [x] **Exercício 2**: Solução numérica para uma fonte de poluição pontual que emite em forma de pulso senoidal $\frac{\partial C}{\partial t}+U \frac{\partial C}{\partial x}=F$ onde $F=F(x,t)$. Os resultados foram apresentados com o diagrama Hovmoller (x, t).
- [x] **Exercício 3**: Similar ao exercício 2 com adição do efeito da difusão: $\frac{\partial C}{\partial t} + U \frac{\partial C}{\partial x} = K \frac{\partial^2C}{\partial x^2} + F$
- [ ] **Exercício 4**: Modelo de água rasa 2D linearizado com fonte de vento zonal tipo ENSO, representativa do efeito de um transiente na intensidade dos ventos alísios. Aplicação da grade tipo C de Arakawa.
- [ ] **Exercício 5**: Aplicação com o modelo de água rasa.



