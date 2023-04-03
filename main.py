import numpy as np
from matplotlib import animation, pyplot as plt
from scipy.spatial import ConvexHull, convex_hull_plot_2d


a, b = 3, 1/2

#q = variable de posición, dq0 = \dot{q}(0) = valor inicial de la derivada
#d = granularidad del parámetro temporal (delta)
def deriv(q,dq0,d):
   #dq = np.empty([len(q)])
   dq = (q[1:len(q)]-q[0:(len(q)-1)])/d
   dq = np.insert(dq,0,dq0) #dq = np.concatenate(([dq0],dq))
   return dq

#Ecuación de un sistema dinámico continuo
#Ejemplo de oscilador simple
def F(q):
    ddq = -8/a*q*(q*q-b)
    return ddq

#Resolución de la ecuación dinámica \ddot{q} = F(q), obteniendo la órbita q(t)
#Los valores iniciales son la posición q0 := q(0) y la derivada dq0 := \dot{q}(0)
def orb(n,q0,dq0,F, args=None, d=0.001):
    q = np.empty([n+1])
    q[0] = q0
    q[1] = q0 + dq0*d
    for i in np.arange(2,n+1):
        args = q[i-2]
        q[i] = - q[i-2] + d**2 * F(args) + 2*q[i-1]
    return q #np.array(q),

def periodos(q,d,max=True):
    #Si max = True, tomamos las ondas a partir de los máximos/picos
    #Si max == False, tomamos las ondas a partir de los mínimos/valles
    epsilon = 5*d
    dq = deriv(q,dq0=None,d=d) #La primera derivada es irrelevante
    if max == True:
        waves = np.where((np.round(dq,int(-np.log10(epsilon))) == 0) & (q >0))
    else:
        waves = np.where((np.round(dq,int(-np.log10(epsilon))) == 0) & (q <0))
    diff_waves = np.diff(waves)
    waves = waves[0][1:][diff_waves[0]>1]
    pers = diff_waves[diff_waves>1]*d
    return pers


#Ejemplo de coordenadas canónicas (q, p)
#Nos quedamos con el más fino y calculamos la coordenada canónica 'p'

def get_pq(q0, dq0, horiz=32):
    # calcular q
    d = 10**(-4)
    n = int(horiz/d)
    t = np.arange(n+1)*d
    q = orb(n,q0=q0,dq0=dq0,F=F,d=d)
    # calcular p
    dq = deriv(q,dq0=dq0,d=d)
    p = dq/2
    return t, p, q

#############################################################################################################
##################                            PROBLEMA 1                               ######################
#############################################################################################################

#Ejemplo de diagrama de fases (q, p) para una órbita completa
def diagrama_de_fases(D0, cmap="winter", horiz=32, linewidth=0.5):
    N = len(D0)
    fig, ax = plt.subplots(1,1, figsize=(12,6))
    colors = [plt.get_cmap(cmap)(i) for i in np.linspace(0,1,N)]
    # plot
    for (q0, p0), col in zip(D0, colors):
        dq0 = 2*p0
        t, p, q = get_pq(q0, dq0, horiz)
        # diagrama de fase de p y q
        ax.plot(q, p, '-', c=col, linewidth=linewidth)
    # settings 
    ax.set_title("Diagrama de fases")
    ax.set_xlabel("q(t)", fontsize=12)
    ax.set_ylabel("p(t)", fontsize=12)
    plt.show()

def problema1():

    x = np.linspace(0., 1., 12)
    y = np.linspace(0., 1., 12)
    D0 = [(x1,x2) for x1 in x for x2 in y]

    diagrama_de_fases(D0, linewidth=2) 

#############################################################################################################
##################                            PROBLEMA 2                               ######################
#############################################################################################################

# A partir de una condiciones iniciales, p0 y q0, calculamos q, p, para cierto t fijo
def get_pq_fijo(t, q0, p0):
    dq0 = 2*p0
    d = 10**(-3)
    n = int(t / d)
    q_ = orb(n,q0=q0,dq0=dq0,F=F,d=d)
    dq_ = deriv(q_,dq0=dq0,d=d)
    q = q_[-1]
    p = dq_[-1] / 2
    return (q,p)

# teniendo en cuenta un cierto D0 -> Calcular D_t para t fijo
def D(D0, t):
    X, Y = [], []
    for q0, p0 in D0:
        q, p = get_pq_fijo(t=t, q0=q0, p0=p0)
        X.append(p)
        Y.append(q)
    return X, Y

# tipo 1
def calcular_area(X, Y):

    ymax, ymin = max(Y), min(Y)
    i1, i2 = Y.index(ymin), Y.index(ymax)
    i1, i2 = min(i1, i2), max(i1, i2)
    d = 10**(-4)
    
    def interseccion(h, Y):
        return np.argmin(abs(Y-h))

    suma = 0
    particion = np.arange(ymin, ymax, step=d)
    error = 0.
    _x1, _x2 = None, None
    for h in particion:
        # dividimos la frontera en dos segmentos: la izq. y der. que van desde 'ymin' hasta 
        # 'ymax', y calculamos la intersección de la recta horizontal con ambos segmentos
        a1 = interseccion(h, np.array(Y[:i1] + Y[i2:])) 
        a2 = interseccion(h, np.array(Y[i1:i2]))
        x1, x2 = X[(a1 if a1 < i1 else a1+(i2-i1))], X[a2+i1]
        suma += abs(x2-x1)
        if _x1 != None: # and _x2 != None
            error += abs(x1-_x1) + abs(x2-_x2)
        _x1, _x2 = x1, x2
    error *= d
    area = suma * d
    return area, error

def plot_conjunto_sumas_riemann(X, Y):

    # separar en dos segmentos
    ymin, ymax = min(Y), max(Y)
    i1, i2 = Y.index(ymin), Y.index(ymax)
    i1, i2 = min(i1, i2), max(i1, i2)
    X1, Y1 = X[:i1] + X[i2:], Y[:i1] + Y[i2:]
    X2, Y2 = X[i1:i2+1], Y[i1:i2+1]

    # plot
    plt.figure(figsize=(5,5))
    plt.plot(X1, Y1, "m-", label="$Y_1$")
    plt.plot(X2, Y2, "g-", label="$Y_2$")
    plt.title("Superficie $D_{t=1/4}$")
    plt.xlabel("$p(t)$")
    plt.ylabel("$q(t)$")
    plt.yticks([ymin, ymax], ["$y_{min}$", "$y_{max}$"])
    plt.xlim(-0.4,1.1)
    plt.ylim(-0.05,1.45)
    plt.legend()
    
    # notas 
    py = Y1[len(Y1)//4]
    px = X1[Y1.index(py)]
    nota = plt.annotate(r'$f_1$',
         xy=(px, py),
         xycoords='data',
         xytext=(0.7, 0.9),
         fontsize=9,
         arrowprops=dict(arrowstyle="->",
         connectionstyle="arc3,rad=.2"))
    py = Y2[len(Y2)//4]
    px = X2[Y2.index(py)]
    nota = plt.annotate(r'$f_2$',
         xy=(px, py),
         xycoords='data',
         xytext=(-0.2, 0.4),
         fontsize=9,
         arrowprops=dict(arrowstyle="->",
         connectionstyle="arc3,rad=.2"))
    
    plt.show()

# tipo 2
def calcular_area_convexhull(seq_q0, seq_dq0, d, title="", plot=False):

    horiz = 0.25

    q2 = np.array([])
    p2 = np.array([])
    for i in range(len(seq_q0)):
        for j in range(len(seq_dq0)):
            q0 = seq_q0[i]
            dq0 = seq_dq0[j]
            n = int(horiz/d)
            q = orb(n,q0=q0,dq0=dq0,F=F,d=d)
            dq = deriv(q,dq0=dq0,d=d)
            p = dq/2
            q2 = np.append(q2,q[-1])
            p2 = np.append(p2,p[-1])

    # convex hull
    X = np.array([q2,p2]).T
    hull = ConvexHull(X)

    if plot: 
        convex_hull_plot_2d(hull)

    return hull.volume

# tipo 2
def calcular_area_2(d, plot=False):
    A_poligonal = calcular_area_convexhull(np.linspace(0.,1.,num=20), np.linspace(0.,2.,num=20), d, title="poligonal", plot=plot)
    A_inferior  = calcular_area_convexhull(np.linspace(0.,1.,num=20), [0], d, title="inferior", plot=plot)
    A_derecha   = calcular_area_convexhull([1], np.linspace(0.,1.,num=20), d, title="derecha")
    A = A_poligonal - A_derecha - A_inferior
    return A
    

def problema2():

    # frontera del conjunto de condiciones iniciales
    border = np.linspace(0.,1.,100) 
    F0 = [(x,0) for x in border] + [(1,x) for x in border] + [(x,1) for x in border[::-1]] + [(0,x) for x in border[::-1]]

    # calculo de D_{1/4}
    X, Y = D(F0, t=0.25)
    
    # area (tipo 1)
    A, error = calcular_area(X, Y)
    plot_conjunto_sumas_riemann(X, Y)
    print("Metodo: Sumas de Riemann")
    print(" - Area: ", A)
    print(" - Error =", error)
    # area (tipo 2)
    Ad4 = calcular_area_2(d=10**(-4), plot=True)
    Ad3 = calcular_area_2(d=10**(-3))
    errord = abs(Ad4 - Ad3)
    print("Metodo: Convex hull")
    print(" - Area: ", Ad4)
    print(" - Error =", errord)

#############################################################################################################
##################                            PROBLEMA 3                               ######################
#############################################################################################################


def f(t, P, Q, ax, colors, markersize, marker):
    p, q = P[:,t], Q[:,t]
    col = colors[t]
    # plot 
    ax.clear()
    e = 0.025
    ax.plot([0,1,1,0,0], [0,0,1,1,0], "-k")
    ax.text(0.5-e, 0.5-e, 
        "$D_0$", fontsize=11, color="black")
    ax.plot(q,p, marker, c=col, markersize=markersize)
    # settings
    e = 0
    ax.set_xlim(-2.5-e,2.5+e)
    ax.set_ylim(-1.5-e,1.5+e)
    ax.set_title("Diagrama de fases para un $t_0$ fijo")
    ax.set_xlabel("$q(t_0)$", fontsize=12)
    ax.set_ylabel("$p(t_0)$", fontsize=12)

def animate_plot(D0, cmap="winter", show=True, namefile="ejemplo.gif", markersize=10, marker="o"):
    P, Q = [], []
    for q0, p0 in D0:
        dq0 = 2*p0
        _, p, q = get_pq(q0, dq0, 5)
        n = len(p)
        k = max(n//100, 1)
        P.append(p[::k])
        Q.append(q[::k])
    N = len(P[0])
    P = np.array(P)
    Q = np.array(Q)
    ts = range(N)
    fig, ax = plt.subplots(1,1, figsize=(12,8))
    colors = [plt.get_cmap(cmap)(i) for i in np.linspace(0,1,N)]
    ani = animation.FuncAnimation(fig, f, ts, fargs=[P, Q, ax, colors, markersize, marker], interval=10)
    ani.save(namefile)
    if show:
        plt.show()

def problema3():

    # frontera F0
    X = np.linspace(0.,1.,100) # 100x4 = 400 puntos
    F0 = [(x,0) for x in X] + [(1,x) for x in X] + [(x,1) for x in X[::-1]] + [(0,x) for x in X[::-1]]

    animate_plot(F0, show=True, namefile="evolucion_diagrama_fases.gif", markersize=5, marker="-")



if __name__ == '__main__':
    while True:
        print("\nElige opcion:")
        print(" 1) Problema 1")
        print(" 2) Problema 2")
        print(" 3) Problema 3")
        print(" 4) Salir")
        x = input(" > ")
        if x == "1":
            problema1()
        elif x == "2":
            problema2()
        elif x == "3":
            problema3()
        else:
            break
