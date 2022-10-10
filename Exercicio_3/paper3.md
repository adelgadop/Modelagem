---
title: "Advecção com difusão e forçante para uma fonte do poluente"
date: \today
author: |
    | Alejandro H. D. Peralta^[Estudante de doutorado, email <aperalta@usp.br>]
    |
    | *Instituto de Astronomia, Geofísica e Ciências Atmosféricas da Universidade de São Paulo*
    | 
keywords: [1D Advecção, Métodos numéricos, RK4]
documentclass: article
link-citations: true
urlcolor: "blue"
classoption: onecolumn
lang: pt
bibliography: "../biblio.bib"
tblPrefix:
  - "Tabela"
  - "Tabelas"
eqnPrefix:
  - "eq."
  - "equações"
csl: "https://raw.githubusercontent.com/citation-style-language/styles/master/aerosol-and-air-quality-research.csl"
abstract: 
    A emissão de um poluente pode variar ao longo do tempo, como no caso de uma chaminé que emite o pulso senoidal no campo básico com velocidade do vento constante. Este trabalho mostra os cálculos para resultados analíticos e numéricos (Euler progressivo-regressivo, leapfrog (2ª e 4ª ordem) e implícito como o esquema Crank-Nicolson). Alguns métodos numéricos geraram oscilações do modo computacional pelo que foram filtrados. Outros experimentos foram considerados para o método implícito com a variação da resolução do tempo $\Delta t$ para obter diferentes números de Courant (CFL) para valores de 1, 2 e 4. A aproximação da ordem 1 é um importante esquema que não precisa de filtros devido à simplicidade do método. No entanto, o esquema é difuso pelo que as concentrações são subestimadas se comparar com a solução analítica. Outros esquemas como leapfrog e Crank-Nicolson geram resultados com oscilações que contradizem o fenômeno físico, pelo que a aplicação de filtros é importante para preservar a monotonicidade. Os resultados dos experimentos são importantes a fim de representar a realidade do fenômeno do transporte dos poluentes na atmosfera, como no caso dos modelos de qualidade do ar.
---

# 1. Introdução
xxx

Conforme com @Doos2020, a discretização de segundo ordem da derivada é expresado como segue,

$$
\left(\frac{d^2 u}{dx^2}\right)_j \approxeq \left[\frac{d}{dx}\left(\frac{d u}{dx}\right) \right]_j \approx \frac{\frac{u_{j+1} - u_j}{\Delta x} - \frac{u_j - u_{j-1}}{\Delta x}}{\Delta x} = \frac{u_{j+1}-2u_j+u_{j-1}}{(\Delta x)^2}
$${#eq:dis}

# 2. Descrição da metodologia
A aproximação considerou as condições do exercício 2 com a adição do efeito da difusão; a equação que governa este problema é dada por:

$$\frac{\partial C}{\partial t} + U \frac{\partial C}{\partial x} = K \frac{\partial^2C}{\partial^2 x} + F$${#eq:ex3}

Onde $F$ é a mesma fonte periódica do Ex. 2, localizado na metade da grade com uma resolução horizontal $\Delta x = 2500$ metros e temporal $\Delta t = 50$ segundos. O requerimento do exercício 3 é determinar o fator $K$ de forma que o tempo de decaimento seja da ordem de 3 horas. Inicialmente $F$ está no tempo n, como segue
$$\frac{\partial C}{\partial t} \rightarrow \frac{(C^{n+1} - C^{n-1})}{2 \Delta t}$$ e radiacional nas condições de fronteira. A @eq:ex3 foi discretizada para o esquema leapfrog (@eq:leap), considerando a advecção e a difusão no tempo n-1 com a forçante no tempo n-1. Depois a forçante é introduzido com o método *splitting*.

$$ C^{n+1}_j =  C^{n-1}_j - 2\Delta t \, U\frac{C^n_{j+1} - C^n_{j-1}}{2\Delta x} + 2\Delta t\,K\frac{C^{n-1}_{j+1}-2C^{n-1}_j+C^{n-1}_{j-1}}{(\Delta x)^2}+ 2\Delta t F^{n-1}_j$${#eq:leap}
ou também expressado como,
$$ C^{n+1}_j =  C^{n-1}_j - \alpha(C^n_{j+1} - C^n_{j-1}) + 2\nu (C^{n-1}_{j+1}-2C^{n-1}_j+C^{n-1}_{j-1})+ 2\Delta t F^{n-1}_j,$$
onde $\alpha = \frac{U\,\Delta t}{\Delta x}$ como número de Courant e $\nu \approxeq K\,\Delta t/(\Delta x)^2$ número de diffusão.

# 3. Resultados


# 4. Discussão dos resultados


# Bibliografia
<div id="refs"></div>

# Apêndice A
xxx



