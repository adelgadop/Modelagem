---
title: "Solução numérica com advecção para o transporte de um poluente "
date: \today
author: "Alejandro H. D. Peralta"
link-citations: true
bibliography: "biblio.bib"

---

# Abstract

O transporte de um poluente pode ser resolvido com diferentes aproximações numéricas. 
Neste trabalho.


# 1. Introdução
A previsão do tempo é importante para entender o impacto da natureza nas atividades humanas e viceversa. Os cientistas desenvolveram soluções baseados na física e matemática. Bjerknes propus sete equações, conhecidas como primitivas, com sete variáveis desconhecidas. Elas requerem de condições de fronteira e iniciais da atmosfera para ser resolvidas. A resolução numérica foi proposta por Richardson com a aplicação de diferenças finitas [@Doos2020] .


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

# 5. Referências