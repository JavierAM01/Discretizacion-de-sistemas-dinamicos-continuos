# Discretización de sistemas dinámicos continuos y teorema de Liouville
	
## Índice

- [Enunciado](#id0)
- [Introducción](#id1)
- [Material usado](#id2)
- [Resultados y conclusiones](#id3)
	- [Pregunta 1](#id3.1)
	- [Pregunta 2](#id3.2)
	- [Pregunta 3](#id3.3)


## Enunciado <a name=id0> </a>

Considera el hamiltoniano de un oscilador no lineal,

$$
H : \mathbb{R}^2 \rightarrow \mathbb{R}, \quad (q,p) \rightarrow p^2 + \dfrac{1}{a}(q^2-b)^2,
$$

donde $a, b\in\mathbb{R}$ son ciertos parámetros. Las ecuaciones de Hamilton–Jacobi del sistema son las
siguientes,

$$
\ddot{q} = -\dfrac{8}{a}q(q^2 - b)
$$

Se trata de una ecuación diferencial que describe la evolución de $q(t)$ y de $p(t) = \dot{q}/2$.

Tomemos $a = 3$ y $b = 1/2$. Supón que disponemos de un conjunto de condiciones iniciales
$D0 := [0, 1] × [0, 1]$, y una granularidad del parámetro temporal $t = n\delta$ con $\delta \in [10−4 , 10−3 ]$,
$n\in \mathbb{N}\cup\{0\}$, con la que puede estimarse la sensibilidad del sistema al grado de discretización.

 1) Representa gráficamente el espacio fásico $D_{(0,\infty)}$ de las órbitas finales del sistema con las
condiciones iniciales $D_0$ . Considera al menos 10 órbitas finales diferentes.

 2) Obtén el valor del área de $D_t$ para $t = 1/4$ y una estimación del su intervalo de error,
presentando los valores de forma cientı́ficamente formal. ¿Se cumple el teorema de Liouville
entre $D_0$ y $D_t$ o bien entre $D_0$ y $D_{(0,\infty)}$?

 3) Realiza una animación GIF con la evolución del diagrama de fases $D_t$ para $t\in (0, 5)$.


## Introducción <a name=id1> </a>
	
La discretización de sistemas dinámicos continuos se refiere al proceso de convertir un sistema dinámico continuo en un sistema discreto, donde el tiempo y las variables se representan en valores discretos. Esta conversión se realiza comúnmente para permitir la simulación y el análisis numérico de sistemas continuos en computadoras digitales. Existen varios métodos de discretización, cada uno tiene sus ventajas y desventajas en términos de precisión, estabilidad y eficiencia.
	
Por otro lado, el teorema de Liouville establece que si se tiene un sistema dinámico Hamiltoniano, entonces la función de densidad de probabilidad en el espacio de fase se conserva a lo largo del tiempo. Esto significa que la probabilidad de que una partícula esté en una región del espacio de fase es constante a medida que el sistema evoluciona.
	
## Material usado <a name=id2> </a>
	
Como lenguaje de programación, se ha usado python, para realizar todo el código, predicciones y gráficas. Para el cálculo del área se ha usado el método de las sumas de Riemann. Como conjunto de datos no se ha usado ningún archivo exterior, únicamente se han utilizado datos creados desde el propio script, com por ejemplo $D_0 = [0,1]\times [0,1]$
	
## Resultados y conclusiones  <a name=id3> </a>
	
### Pregunta 1  <a name=id3.1> </a>
	
Escogemos una serie de puntos aleatorios en el conjunto de condiciones iniciales, $(q(0), p(0)) \in D_0 := [0,1]\times [0,1]$. Para cada par $(q(0), p(0))$ calculamos y graficamos el espacio fásico $D_{(0,\infty)}$ de la órbitas finales del sistema. 


<div style="text-align:center;">
  <image src="/images/p1_1.png" style="width:70%; height:8cm;" alt="Espacio fásico de la órbitas finales">
</div>

### Pregunta 2  <a name=id3.2> </a>

Para obtener el valor del área se ha procesado mediante sumas de Riemann. Para ello se ha recorrido desde el mínimo valor en $q$, $y_{min}$, hasta el máximo, $y_{max}$. Generamos una partición de dicho intervalo, $\mathcal{P} = \{ h_0 = y_{min}, h_1, \dots, h_n=y_{max} \}$. Donde $h_{i+1} - h_i = d = 10^{-4}$ es fijo para todo $i = 0,\dots, n-1$. Para calcular el error cometido, calculamos el error por cada intervalo, $e_i = d\cdot (x_{i+1}-x_i)$ y así el error final es, $\displaystyle{E = \sum_{i=1}^n e_i = 0.000246}$. Finalmente el área final es: $A = 1.0006 \pm (2)$.

Hay que tener en cuenta que la función a integrar no es continua, es discreta. Por tanto no se puede calcular directamente la suma superior o inferior de las sumas de Rieman. En este caso se está integrando desde el eje $y$ (variable $q$) ya que las dos funciones $f_1, f_2:[y_{min}, y_{max}] \rightarrow \mathbb{R}$ son inyectivas. Para hallar la altura del rectángulo del subintervalo $[h_i,h_{i+1}] \subset [y_{min}, y_{max}]$ se obtiene el valor más cercano a $h_{i+1}$, $y^{i+1}$, y luego se obtiene su altura respectiva, $x_j^{i+1}$ para $j=1,2.$, es decir, $x_j^{i+1} = X[a_j^{i+1}]$ donde $a_j^{i+1} = \underset{y\in Y_j}{\arg \min} [\text{abs}(y-h_{i+1})]$, siendo $Y_j = f_j([y_{min}, y_{max}])$, para $j=1,2$.

<div style="text-align:center;">
  <image src="/images/p2_1.jpg" style="width:45%; height:8cm;" alt="Sumas de Riemann">
  <image src="/images/p2_2.png" style="width:35%; height:8cm;" alt="Conjunto de condiciones inciales, para $t=\frac{1}{4}$">
</div>

Por el teorema de Liouville, sabemos que el área comprendido por el espacio fásico es constante. Además $D_0 := [0,1]\times [0,1]$ por tanto es fácil ver que su área es $A = 1$. Por lo tanto, por el teorema de Liouville, el área de $D_t$ es igual a 1, para cualquier $t\in \mathbb{R^+}$. Efectivamente, para $t=\frac{1}{4}$ nos ha dado $A = 1.0006 \pm (2)$ por lo que es una buena aproximación del valor real, cumpliendose así el teorema de Liouville.

También podríamos haber usado la función de convexhull, dada en la plantilla. Nuevamente nos da valores muy cercanos a 1, en concreto se ha calculado el área para $d=10^{-4}$ y $d=10^{-3}$ y restando ambas se ha caculado el error, dando así un área de: $A = 1.0044 \pm (5)$. Hay que tener en cuenta que el conjunto no es convexo por lo que habría que restar fragmentos del conjunto. En todo caso se cumple el Teorema de Liouville.

<div style="text-align:center;">
  <image src="/images/ch.png" style="width:40%; height:8cm;" alt="Conjunto $D_{t=1/4}$ no convexo">
  <image src="/images/ch2.png" style="width:40%; height:8cm;" alt="Ejemplo de fragmento a restar">
</div>
 
### Pregunta 3  <a name=id3.3> </a>

Para realizar una animación gif con la evolución del diagrama de fases para $t\in (0,5)$, procemos como el ejercicio anterior para calcular $D_{t=1/4}$ pero graficando múltiples $t\in (0,5)$. Para realizar la figura se ha graficado únicamente su frontera. Se pueden graficar también los puntos del interior pero la frontera se aprecia menos. Aunque también, al graficar únicamente la frontera con un conjunto finito de puntos se pueden apreciar imprecisiones en la unión de los mismos en el núcleo de la espiral. El gif se deja adjunto como *evolucion_diagrama_fases.gif*.

<div style="text-align:center;">
  <img src="/images/evolucion_diagrama_fases.gif" style="width:100%; height:12cm;">
</div>

