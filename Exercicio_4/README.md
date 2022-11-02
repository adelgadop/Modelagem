## Exercício 4 - Água Rasa 2D
Entrega do trabalho  + apresentações nos dias 28/10/2022

---

![Agua Rasa, H=50 m, plano B](Exercicio_4/gifs/leap_50m_scen3.gif)

O modelo de água rasa é usado na modelagem atmosférica e oceanográfica para uma enorme variedade de problemas.  Na atmosfera, por ex., uma fonte de calor na atmosfera (por exemplo, um furacão) pode ser vista como uma fonte de massa no modelo de água rasa.  No oceano, o gatilho para formação do ENSO pode ser representado por uma transferência de momentum na componente zonal do vento da região equatorial, representativa do efeito de um transiente na intensidade dos ventos alísios.  Como esses problemas são altamente dependentes da capacidade do modelo descrever o processo de dispersão de energia, vamos fazer simulações com o modelo de água rasa 2D:

$\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y} - fv + \frac{\partial \phi}{\partial x} = F_u$

$\frac{\partial v}{\partial t}+u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} - fu + \frac{\partial \phi}{\partial y} = F_v$

$\frac{\partial \phi}{\partial t}+u\frac{\partial \phi}{\partial x} + v\frac{\partial \phi}{\partial y} +\phi.Div(V) = F_{\phi}$

Onde $F_u$, $F_v$, $F_{\phi}$ representam fontes de momentum zonal, meridional e massa, respectivamente.

a. Resolva em 2D o modelo de água rasa (1) na versão linear, com estado básico U = 0, supondo que 
- i. f = 0
- f = f$_o$ = constante (assuma que a latitude é 20°S)
- f = $\beta$y (ou seja, o plano beta equatorial). Suponha que a latitude máxima seja de +/- 4000 km, com dx=dy=100 km. Na longitude também considere o domínio +/- 4000 km.
- Opções de fonte (um problema meteorológico e outro oceanográfico) - escolha uma deles:
   - $\bf A.$ Aplique a fonte de massa na atmosfera constante no tempo. A estrutura espacial da fonte é dada por uma Gaussiana (use Nr$_x$=Nr$_y$= 4 e 10 - isto é, fonte pequena e fonte grande). Considere que a fonte esteja centrada na latitude de +1500 km e longitude 0 (centro da malha em x). Use o domínio horizontal de 8000 km. Use CF aberta radiacional nas fronteiras. Utilize `H=250` m e integre por 3 dias.
   - $\bf B.$ Fonte de momentum zonal, constante (vento forçante de leste), com gaussiana centrada no equador, alongada na direção zonal. (Nr$_x$=10 e Nr$_y$=4) para representar a troca de momentum (representado o processo de formação do ENSO). Centre a fonte de momentum na longitude 0. Use CF radiacional na fronteira oeste, norte e sul e CF rígida na fronteira leste. Use `H = 1 m` (sim, um metro!). Integre até que a onda que vai para leste atinja a costa e comece a refletir para o centro da grade. Como vai demorar muito tempo, faça os testes do modelo com `H = 250 m`. Só depois que você tiver removido todas as "funcionalidades não documentadas" do modelo (i.e., os "bugs"), faça uma integração com o H mais realista.

b. Discuta os resultados (do ponto de vista da dispersão de energia) e as diferenças entre os métodos numéricos. Use pelo menos 2 esquemas numéricos de integração no tempo (um deles deve ser o leapfrog). Use a grade C.

c. Em todos os casos acima, verifique a conservação de massa e de energia e faça mapas da divergência e vorticidade. Discuta as diferenças entre os métodos numéricos.

d. Observe que como o modelo é linear, a magnitude da fonte de massa ou de momentum (dada pela amplitude da gaussiana) é irrelevante. Considere a fonte com magnitude máxima unitária.

e. E divirta-se...

Faça o relatório do Ex. 4 na forma de artigo.






