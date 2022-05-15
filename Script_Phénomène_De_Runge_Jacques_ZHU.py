#Auteur : Jacques ZHU

    #Modules
import numpy as np
import matplotlib.pyplot as plt

    #Variables
n=10   #nombre de points xi
N= 100  #résolution (dx)
a=-1; b=-a  #bornes de l'étude fonctionnelle
x= 0.5 #abscisse


    #Fonctions set de points

def x1(i,N):   #Points -1+2i/N avec N le nombre de points
    return -1+2*i/N

def x2(i,N):   #Points cos((2i+1)pi/(2N+2)) avec N le nombre de points
    return np.cos((2*i+1)*np.pi/(2*N+2))


def Liste_xi(N,fxi):  #Liste des N+1 points xi entre 0 et N selon la fonction set de points fxi
    L=[]
    for i in range (0,N+1):
        L.append(fxi(i,N))
    return L
#print(Liste_xi(10,x1))


    #Fonctions à interpoler

def g(x):
    return 1/(1+25*x*x)

def identite(x):
    return 1


    #Interpolation de Lagrange

def monome(x,x0): #Polynome de degré 1 qui s'annule en x0 de coeff directeur 1
    return x-x0

def Li(x,xi,Lxi):   #Polynome interpolateur élémentaire i de Lagrange évalué en x et en la liste Lxi des points
    P=1
    for j in range (0,len(Lxi)):
        if xi==Lxi[j]:
            P=P
        else:
            P*=monome(x,Lxi[j])/monome(xi,Lxi[j])
    return P


def P_Lagrange(f,fxi,x,n): #Polynome d'interpolation de Lagrange de f interpolé par le set de points de fonction fxi à n points en x
    Lxi=Liste_xi(n,fxi)
    y=0
    for i in range (0,len(Lxi)):
        y+=f(Lxi[i])*Li(x,Lxi[i],Lxi)
    return y

#print(P_Lagrange(g,x1,x,n))
#print(P_Lagrange(g,x2,x,n))


def P_LagrangeAbs(f,fxi,x,n): #Polynome d'interpolation de Lagrange avec la valeur absolue des élémentaires
    Lxi=Liste_xi(n,fxi)
    y=0
    for i in range (0,len(Lxi)):
        y+=f(Lxi[i])*abs(Li(x,Lxi[i],Lxi))
    return y

#print(P_LagrangeAbs(g,n,x1,x))
#print(P_LagrangeAbs(g,n,x2,x))


    #Graph de l'interpolation de Lagrange

def liste(f,a,b,fxi,n,N): #Création de 3 listes de taille N des points d'abscisses, d'images de la fonction et d'images de son interpolation entre a et b de degré n
    Lx = np.linspace(a,b,N)
    Lf1= [f(i) for i in Lx]
    Lf2= [P_Lagrange(f,fxi,x,n) for x in Lx]
    return Lx,Lf1,Lf2

def graph(fu,a,b,fxi,n,N):  #Courbes de la fonction et de l'interpolation entre a et b avec n points de précision N
    Lx, Lf1, Lf2 = liste(fu,a,b,fxi,n,N)
    plt.figure()
    plt.plot(Lx, Lf1)
    plt.plot(Lx, Lf2)
    return plt.show()
#graph(g,a,b,x1,10,100)
#graph(g,a,b,x2,15,200)


def graphN(fu,a,b,fxi,n0,n,N):  #Courbes de la fonction et de l'interpolation entre a et b de n0 jusqu'à n points
    plt.figure()
    Lx = np.linspace(a,b,N)
    Lf1= [fu(i) for i in Lx]
    plt.plot(Lx, Lf1,label="g")
    for i in range (n0,n+1):
        Lf2= [P_Lagrange(fu,fxi,x,i) for x in Lx]
        C='Polynome n= '+str(i)
        plt.plot(Lx,Lf2,label =C)
    plt.legend()
    return plt.show()

#graphN(g,a,b,x1,7,10,100)
#graphN(g,a,b,x2,7,10,100)

def graphN2(fu,a,b,fxi,L,N):  #Courbes de la fonction et de l'interpolation entre a et b avec n points avec n dans L
    plt.figure()
    Lx = np.linspace(a,b,N)
    Lf1= [fu(i) for i in Lx]
    plt.plot(Lx, Lf1,label="g")
    for i in L:
        Lf2= [P_Lagrange(fu,fxi,x,i) for x in Lx]
        C='Polynome n= '+str(i)
        plt.plot(Lx,Lf2,label =C)
    plt.legend()
    return plt.show()
#graphN2(g,a,b,x1,[5,7,9],500)
#graphN2(g,a,b,x2,[7,10,20,25],500)


    #Graphe d'erreur entre f et Pn (Norme infinie de la différence entre f et Pn)

def sup(L1,L2):   #Obtient le maximum d'écart d'éléments entre 2 listes
    M=0
    for i in range (len(L1)):
        M = max([M,abs(L1[i]-L2[i])])
    return M

def supn(fu,a,b,n,precision,points):  #Obtient une courbe du maximum d'écart entre la fonction et son interpolation (x in [a,b]) en fonction du nombre de points interpolés
    Ln= []
    SupDf = []
    for i in range (1,n):
        Ln.append(i)
        inutile,Lf1,Lf2=liste(fu,a,b,points,i,precision)
        SupDf.append(sup(Lf1,Lf2))
    plt.figure()
    plt.plot(Ln,SupDf, label='||g-Pn||')
    plt.grid()
    plt.show()
    return 0

#supn(g,a,b,10,100,x1)
#supn(g,a,b,10,100,x2)


#Calcul de Lambda(n)
def lambd(n,x,points): #calcule la petite constante de Lambda dépendant de n en x
    return P_LagrangeAbs(identite,points,x,n), P_Lagrange(identite,points,x,n)


def supLambda(a,b,n,precision,points): #calcule le sup en fonction de n : renvoie le graphe Lambda(n) qui dépend des points
    Ln = [i for i in range (1,n+1)]
    dx = np.linspace(a,b,precision)
    Lmax = []
    for i in range(1,n+1) :
        max=0
        for x in dx:
            lamb=lambd(i,x,points)[0]
            if max<lamb:
                max= lamb
        Lmax.append(max)
    plt.figure()
    plt.plot(Ln,Lmax, label='Lambda(n,famille_xi)')
    plt.legend()
    return plt.show()

#supLambda(-1,1,15,500,x1)
#supLambda(-1,1,15,500,x2)





#Figures

    #Figure 1 Introduction : Queue de Chat

# Taille = [511,287]
# ListeX = [131,145,160,166,174,182,189,206,222,239,255,275,289,298,308,322,340,356,368]
# ListeYBrut = [101,117,129,133,137,141,145,151,158,161,164,166,169,172,175,178,185,192,199]
# ListeY = [Taille[1]-y for y in ListeYBrut]
# def P_Lagrange_Points(Lx,Ly,x): #Polynome d'interpolation de Lagrange de f interpolé par fxi a n points en x
#     y=0
#     for i in range (0,len(Lx)):
#         y+=Ly[i]*Li(x,Lx[i],Lx)
#     return y
# Queue_De_Chat = (np.linspace(0,511,512), [P_Lagrange_Points(ListeX,ListeY,x) for x in np.linspace(0,511,512)])
# plt.figure()
# plt.plot(Queue_De_Chat[0],Queue_De_Chat[1])
# plt.show()


    #Figure 2

#graphN(g,a,b,x1,1,6,100)
#graphN(g,a,b,x1,7,10,100)


    #Figure 3

#supn(g,a,b,20,200,x1)


    #set de Points équirépartis et Tchebychev

# Lx1 = Liste_xi(50,x1)
# Lpas1= [0 for x in Lx1]
# Lx2 = Liste_xi(50,x2)
# Lpas2= [-1 for x in Lx2]
# plt.figure()
# plt.plot(Lx1,Lpas1,label='points équirépartis', marker='o',markersize=1, ls='None')
# plt.plot(Lx2,Lpas2,label='points de Tchebychev',marker='o',markersize=1, ls='None')
# plt.legend()
# plt.show()


    #Figure 4

#supLambda(-1,1,15,500,x1)
#supLambda(-1,1,15,500,x2)


    #Figure 5
#graphN2(g,a,b,x2,[7,10,20,25],500)
#supn(g,a,b,30,500,x2)