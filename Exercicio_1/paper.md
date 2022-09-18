---
title: "Solução numérica com advecção para o transporte de um poluente "
date: \today
author: "Alejandro H. D. Peralta"
link-citations: true
urlcolor: "blue"
bibliography: "biblio.bib"
csl: "https://raw.githubusercontent.com/citation-style-language/styles/master/aerosol-and-air-quality-research.csl"
---

# Abstract

O transporte de um poluente pode ser resolvido com diferentes aproximações numéricas. 
Neste trabalho.


# 1.Introdução
A previsão do tempo é importante para entender o impacto da natureza nas atividades humanas e viceversa. Os cientistas desenvolveram soluções baseados na física e matemática. Bjerknes propus sete equações, conhecidas como primitivas, com sete variáveis desconhecidas [@Alvim2013]. Elas requerem de condições de fronteira e iniciais da atmosfera para ser resolvidas. A resolução numérica foi proposta por Richardson com a aplicação de diferenças finitas [@Doos2020] .

A condição inicial (CI) dada por uma gaussiana centrada em i=51, com decaimiento exponencial dado por `nr` (número de pontos) onde a amplitude da perturbação cai de um fator e:

$$
C(x,0) = C_{i,0} = C_0 \exp[\frac{-(i\Delta x - 51\Delta x)^2}{(nr*\Delta x)^2}]
$$

CFL:
$$
U*\frac{\Delta t}{\Delta x} < 1
$$
## 1.1 Aproximação Leapfrog
Discretização no tempo e espaço das derivadas considerando diferenças centrais como segue:

$$
\frac{C^{n+1}_i - C^{n-1}_i}{2\Delta t} = -U*\frac{C^n_{i+1} - C^n_{i-1}}{2\Delta x}
$$

Se CFL é igual a $\gamma$, então temos

$$
C^{n+1}_i = C^{n-1}_i - \gamma(C^n_{i+1} - C^n_{i-1})
$$

para i = 0, consideramos usar o esquema de aproximação progressiva central no tempo:

$$
\frac{C^{n+1}_i - C^{n}_i}{\Delta t} = -U*\frac{C^n_{i+1} - C^n_{i-1}}{2\Delta x}
$$

Com isso temos

$$
C^{1}_i = C^{0}_i - \frac{\gamma}{2}(C^0_{i+1} - C^0_{i-1})
$$

# 2. Descrição da metodologia
Temos uma equação da adveção em 1D:

$$
\frac{\partial C}{\partial t} + U*\frac{\partial C}{\partial x} = 0,
$$

Temos três maneiras possíveis para representar a forma discreta de 
$$
\frac{\partial C}{\partial x},
$$

- Diferenças progressivas
- Diferenças regressivas
- Diferenças centradas

O esquema de diferenças finitas escolhido "progressivo no tempo e regressivo no espaço" é um método de ordem 1.

# 3. Resultados
xxxxx

# 4. Discussão dos resultados

# Bibliografia