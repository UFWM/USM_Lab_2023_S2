import matplotlib.pyplot as plt 
import numpy as np 
import math 
import random

def ALG(D,Q,p,C,r):

  if(Q[0]<D[0] or Q[1]+(Q[0]-D[0])*(1-r)<D[1] or
     Q[2]+(Q[1]-D[1])*(1-r)+(Q[0]-D[0]-((D[1]-Q[1])/(1-r)))*(1-r)<D[2]):
    return [C[2][2]*Q[2]+C[1][2]*Q[1]+C[0][2]*Q[0],0,0,0,0,0,0,0,0,0]
  
  eta=0
  if D[1]>Q[1]:
    eta=1

  prodMax=[0.0,0.0,0.0]
  prod=[0.0,0.0,0.0]
  flujo=[0.0,0.0,0.0]

  delta=[0.0,0.0,0.0]
  DELTA=[0.0,0.0,0.0]
  limProd=[[Q[0]*p[0][0],0.0,0.0],[Q[1]*p[1][0],0.0,0.0],[Q[2]*p[2][0],0.0,0.0]]
  for i in range(3):
    for j in range(2):
      limProd[i][j+1]=limProd[i][j]+p[i][j+1]*Q[i]

  for i in range(3):
    prod[i]=min(D[i],limProd[i][0])
    prodMax[i]=prod[i]

  if(prodMax[1]>=D[1]):
    DELTA[1]=limProd[1][0]-D[1]
    delta[2]=D[2]-prodMax[2]
    if(prodMax[2]<D[2] and prodMax[2]+DELTA[1]*(1-r)>=D[2]):
      flujo[1]+=delta[2]/(1-r)
      prod[1]+=delta[2]/(1-r)
      prodMax[2]=D[2]
    if(prodMax[2]<D[2] and prodMax[2]+DELTA[1]*(1-r)<D[2]):
      flujo[1]+=DELTA[1]
      prod[1]+=DELTA[1]
      prodMax[2]+=DELTA[1]*(1-r)

  if(prodMax[0]>=D[0]):
    DELTA[0]=limProd[0][0]-prod[0]
    delta[2]=D[2]-prodMax[2]
    if(limProd[2]<D[2] and limProd[2]+DELTA[0]*(1-r)>=D[2] and (Q[0]-prod[0]-delta[2]/(1-r))*(1-r)>=D[1]):
      flujo[2]+=delta[2]/(1-r)
      prod[0]+=delta[2]/(1-r)
      prodMax[2]+=D[2]
      DELTA[0]-=delta[2]/(1-r)
    delta[1]=D[1]-prodMax[1]
    if(prodMax[1]<D[1] and prodMax[1]+DELTA[0]*(1-r)>=D[1]):
      flujo[0]+=delta[1]/(1-r)
      prod[0]+=delta[1]/(1-r)
      prodMax[1]+=D[1]
      DELTA[0]-=delta[1]/(1-r)
    if(limProd[2]<D[2] and limProd[2]+DELTA[0]*(1-r)<D[2] and (Q[0]-prod[0]-delta[2]/(1-r))*(1-r)>=D[1]):
      flujo[2]+=DELTA[0]
      prod[0]+=DELTA[0]
      prodMax[2]+=DELTA[0]*(1-r)
      DELTA[0]=0
    if(prodMax[1]<D[1] and prodMax[1]+DELTA[0]*(1-r)<D[1]):
      flujo[0]+=DELTA[0]
      prod[0]+=DELTA[0]
      prodMax[1]+=DELTA[0]*(1-r)

  if(prodMax[1]<D[1]):
    m=min(limProd[1][1]-limProd[1][0],D[1]-prodMax[1])
    prod[1]+=m
    prodMax[1]+=m

  if(prodMax[0]<D[0]):
    m=min(limProd[0][1],D[0])
    prod[0]=m
    prodMax[0]=m

  if(prodMax[2]<D[2]):
    m=min(limProd[2][1]-limProd[2][0],D[2]-prodMax[2])
    prod[2]+=m
    prodMax[2]+=m

  if(prodMax[1]>=D[1]):
    DELTA[1]=limProd[1][1]-prod[1]
    delta[2]=D[2]-prodMax[2]
    if(prodMax[2]<D[2] and prodMax[2]+DELTA[1]*(1-r)>=D[2]):
      flujo[1]=flujo[1]+delta[2]/(1-r)
      prod[1]=prod[1]+delta[2]/(1-r)
      prodMax[2]=D[2]
    if(prodMax[2]<D[2] and prodMax[2]+DELTA[1]*(1-r)<D[2]):
      flujo[1]=flujo[1]+DELTA[1]
      prod[1]=prod[1]+DELTA[1]
      prodMax[2]=prodMax[2]+DELTA[1]*(1-r)

  if(prodMax[0]>=D[0]):
    DELTA[0]=limProd[0][1]-prod[0]
    delta[2]=D[2]-prodMax[2]
    if(prodMax[2]<D[2] and prodMax[2]+DELTA[0]*(1-r)>=D[2] and C[0][1]<=C[2][2]*(1-r) and (Q[0]-prod[0]-delta[2]/(1-r))*(1-r)>=D[1]):
      flujo[2]+=delta[2]/(1-r)
      prod[0]+=delta[2]/(1-r)
      prodMax[2]+=D[2]
      DELTA[0]-=delta[2]/(1-r)
    delta[1]=D[1]-prodMax[1]
    if(prodMax[1]<D[1] and C[0][1]<=C[2][2]*(1-r) and prodMax[1]+DELTA[0]*(1-r)>=D[1]):
      flujo[0]+=delta[1]/(1-r)
      prod[0]+=delta[1]/(1-r)
      prodMax[1]+=D[1]
      DELTA[0]-=delta[1]/(1-r)
    if(prodMax[2]<D[2] and prodMax[2]+DELTA[0]*(1-r)<D[2] and (Q[0]-prod[0]-delta[2]/(1-r))*(1-r)>=D[1]):
      flujo[2]+=DELTA[0]
      prod[0]+=DELTA[0]
      prodMax[2]+=DELTA[0]*(1-r)
      DELTA[0]=0
    if(prodMax[1]<D[1] and prodMax[1]+DELTA[0]*(1-r)<D[1]):
      flujo[0]+=DELTA[0]
      prod[0]+=DELTA[0]
      prodMax[1]+=DELTA[0]*(1-r)

  if(prodMax[2]<D[2]):
    m=min(limProd[2][2]-limProd[2][1],D[2]-prodMax[2])
    prod[2]=prod[2]+m
    prodMax[2]=prodMax[2]+m

  if(prodMax[1]<D[1]):
    m=min(limProd[1][2]-limProd[1][1],D[1]-prodMax[1])
    prod[1]=prod[1]+m
    prodMax[1]=prodMax[1]+m

  if(prodMax[0]<D[0]):
    prod[0]+=D[0]-prodMax[0]
    prodMax[0]=D[0]

  if(prodMax[1]>=D[1] and prodMax[2]<D[2]):
    DELTA[1]=limProd[1][2]-prod[1]
    if(prodMax[2]<D[2] and prodMax[2]+DELTA[1]*(1-r)>=D[2]):
      delta[2]=D[2]-limProd[2]
      flujo[1]+=delta[2]/(1-r)
      prod[1]+=delta[2]/(1-r)
      prodMax[2]=D[2]
    if(prodMax[2]<D[2] and prodMax[2]+DELTA[1]*(1-r)<D[2]):
      delta[2]=D[2]-limProd[2]
      flujo[1]+=DELTA[1]
      prod[1]+=DELTA[1]
      prodMax[2]=prodMax[2]+DELTA[1]*(1-r)

  if(prodMax[0]>=D[0] and prodMax[2]<D[2]):
    DELTA[0]=limProd[0][2]-prod[0]
    if(prodMax[2]+DELTA[0]*(1-r)>=D[2]):
      delta[2]=D[2]-prodMax[2]
      flujo[2]+=delta[2]/(1-r)
      prod[0]+=delta[2]/(1-r)
      prodMax[2]=D[2]
    if(prodMax[2]+DELTA[0]*(1-r)<D[2]):
      flujo[2]+=DELTA[0]
      prod[0]+=DELTA[0]
      prodMax[2]+=DELTA[0]*(1-r)

  if(prodMax[0]>=D[0] and prodMax[1]<D[1]):
    DELTA[0]=limProd[0][2]-prod[0]
    if(prodMax[1]+DELTA[0]*(1-r)>=D[1]):
      delta[1]=D[1]-prodMax[1]
      flujo[0]+=delta[1]/(1-r)
      prod[0]+=delta[1]/(1-r)
      prodMax[1]=D[1]
    if(prodMax[1]+DELTA[0]*(1-r)<D[1]):
      flujo[0]+=DELTA[0]
      prod[0]+=DELTA[0]
      prodMax[1]+=DELTA[0]*(1-r)

  V=0
  for i in range(3):
    if(prod[i]>limProd[i][0] and prod[i]<=limProd[i][1]):
      V=V+C[i][1]*(prod[i]-limProd[i][0])
    if(prod[i]>limProd[i][1]):
      V=V+C[i][2]*(prod[i]-limProd[i][1])+C[i][1]*(limProd[i][1]-limProd[i][0])

  resAlg=[round(V,2),round(prod[0],2),round(prod[1],2),round(prod[2],2),round(flujo[0],2),0.0,round(flujo[1],2),0.0,round(flujo[2],2),0.0]
  return resAlg
def permutarALG(D,Q,p,C,r):

  perm=[]
  perm.append(ALG([D[0],D[1],D[2]],[Q[0],Q[1],Q[2]],[p[0],p[1],p[2]],[C[0],C[1],C[2]],r))
  perm.append(ALG([D[0],D[2],D[1]],[Q[0],Q[2],Q[1]],[p[0],p[2],p[1]],[C[0],C[2],C[1]],r))
  perm.append(ALG([D[1],D[0],D[2]],[Q[1],Q[0],Q[2]],[p[1],p[0],p[2]],[C[1],C[0],C[2]],r))
  perm.append(ALG([D[1],D[2],D[0]],[Q[1],Q[2],Q[0]],[p[1],p[2],p[0]],[C[1],C[2],C[0]],r))
  perm.append(ALG([D[2],D[0],D[1]],[Q[2],Q[0],Q[1]],[p[2],p[0],p[1]],[C[2],C[0],C[1]],r))
  perm.append(ALG([D[2],D[1],D[0]],[Q[2],Q[1],Q[0]],[p[2],p[1],p[0]],[C[2],C[1],C[0]],r))

  resultadoPerm=perm[0]
  permSol=1
  for j in range(5):
    if resultadoPerm[0]>perm[j+1][0]:
      resultadoPerm=perm[j+1]
      permSol=j+2

  #print("\n\nPERM: ",permSol,", Demandas: ",D[0],D[1],D[2],"\n")
  #printSol(resultadoPerm)
  ##RESULTADO(valorObj, Norte, Sur, Centro, NS, SN, SC, CS, NC, CN)
  #print(permSol)

  if permSol==1:
    ##NSC
    #print("PERM 1:\nD: ",D[0],D[1],D[2],"\nQ: ",Q[0],Q[1],Q[2],"\nC: ",C[0],C[1],C[2],"\np: ",p[0],p[1],p[2],"\nr: ",r)
    return resultadoPerm
  elif permSol==2:
    ##NCS
    #print("PERM 2:\nD: ",D[0],D[2],D[1],"\nQ: ",Q[0],Q[2],Q[1],"\nC: ",C[0],C[2],C[1],"\np: ",p[0],p[2],p[1],"\nr: ",r)
    return [resultadoPerm[0],resultadoPerm[1],resultadoPerm[3],resultadoPerm[2],
            resultadoPerm[8],resultadoPerm[9],resultadoPerm[7],resultadoPerm[6],resultadoPerm[4],resultadoPerm[5]]
  elif permSol==3:
    ##SNC
    #print("PERM 3:\nD: ",D[1],D[0],D[2],"\nQ: ",Q[1],Q[0],Q[2],"\nC: ",C[1],C[0],C[2],"\np: ",p[1],p[0],p[2],"\nr: ",r)
    return [resultadoPerm[0],resultadoPerm[2],resultadoPerm[1],resultadoPerm[3],
            resultadoPerm[5],resultadoPerm[4],resultadoPerm[8],resultadoPerm[9],resultadoPerm[6],resultadoPerm[7]]
  elif permSol==4:
    ##SCN
    #print("PERM 4:\nD: ",D[1],D[2],D[0],"\nQ: ",Q[1],Q[2],Q[0],"\nC: ",C[1],C[2],C[0],"\np: ",p[1],p[2],p[0],"\nr: ",r)
    return [resultadoPerm[0],resultadoPerm[2],resultadoPerm[3],resultadoPerm[1],
            resultadoPerm[6],resultadoPerm[7],resultadoPerm[9],resultadoPerm[8],resultadoPerm[5],resultadoPerm[4]]
  elif permSol==5:
    ##CNS
    #print("PERM 5:\nD: ",D[2],D[0],D[1],"\nD: ",D[2],D[0],D[1],"\nC: ",C[2],C[0],C[1],"\np: ",p[2],p[0],p[1],"\nr: ",r)
    return [resultadoPerm[0],resultadoPerm[3],resultadoPerm[1],resultadoPerm[2],
            resultadoPerm[9],resultadoPerm[8],resultadoPerm[4],resultadoPerm[5],resultadoPerm[7],resultadoPerm[6]]
  elif permSol==6:
    ##CSN
    #print("PERM 6:\nD: ",D[2],D[1],D[0],"\nD: ",D[2],D[1],D[0],"\nC: ",C[2],C[1],C[0],"\np: ",p[2],p[1],p[0],"\nr: ",r)
    return [resultadoPerm[0],resultadoPerm[3],resultadoPerm[2],resultadoPerm[1],
            resultadoPerm[7],resultadoPerm[6],resultadoPerm[5],resultadoPerm[4],resultadoPerm[9],resultadoPerm[8]]
  
  return resultadoPerm
def interpol(interpolHoraN,interpolHoraS,interpolHoraC,n,var):
  D=[]
  for i in range(len(interpolHoraN)-1):
    for j in range(n):
      demInstN=interpolHoraN[i]+(j/n)*(interpolHoraN[i+1]-interpolHoraN[i])
      demInstS=interpolHoraS[i]+(j/n)*(interpolHoraS[i+1]-interpolHoraS[i])
      demInstC=interpolHoraC[i]+(j/n)*(interpolHoraC[i+1]-interpolHoraC[i])
      D.append([max([round(demInstN*(1+var*(2*random.random()-1))),0]),max([round(demInstS*(1+var*(2*random.random()-1))),0]),max([round(demInstC*(1+var*(2*random.random()-1))),0])])
  return D
def printSol(respuesta):
  print("=====Resultados=====")
  print("")
  print("Costo minimo: $", respuesta[0])
  print("Produccion Norte: ", respuesta[1])
  print("Produccion Sur: ", respuesta[2])
  print("Produccion Centro: ", respuesta[3])
  print("")
  print("Flujo NS: ", respuesta[4])
  print("Flujo SN: ", respuesta[5])
  print("Flujo SC: ", respuesta[6])
  print("Flujo CS: ", respuesta[7])
  print("Flujo NC: ", respuesta[8])
  print("Flujo CN: ", respuesta[9])
  print("")

#=====Parametros=====

random.seed(0) #Quita aleatoridad
r=0.005 #0.5%
#Interpolaciones de los nodos N,S,C
interpolHoraN=[10.60,10.00,9.62,9.60,9.80,11.00,13.54,14.80,15.20,15.41,15.50,15.45,15.20,15.10,15.00,14.80,14.90,15.35,16.04,15.90,15.50,14.60,13.10,11.80]
interpolHoraS=[180,160,150,160,220,260,370,500,580,560,530,540,530,510,500,520,550,620,640,630,500,420,240,200]
interpolHoraC=[0.4,0.37,0.34,0.33,0.35,0.43,0.55,0.74,0.85,0.82,0.80,0.77,0.78,0.75,0.74,0.75,0.78,0.90,0.92,0.86,0.73,0.60,0.50,0.43]
#Demanda promedio esperada
demPromedioTotal=10500
#Factor de aporte de cada nodo al promedio total
pondNSC=[0.3,0.1,0.6]
#n-1 : Puntos a interpolar linealmente entre cada dato de la interpolacion
n=12
#var : Porcentaje de variacion de los datos aplicada a la interpolacion
var=0.02
#====================

#Arreglos
resultado, listaD, Q, p, C=[], [], [6000,2000,12000], [[0.3,0.3,0.4],[0.1,0.4,0.5],[0.2,0.1,0.7]], [[0,40,80],[0,40,80],[0,40,80]]

#Calculo interpolacion
promInterpolN=sum(interpolHoraN)/len(interpolHoraN)
promInterpolS=sum(interpolHoraS)/len(interpolHoraS)
promInterpolC=sum(interpolHoraC)/len(interpolHoraC)

for i in range(len(interpolHoraN)):
  interpolHoraN[i]=round((interpolHoraN[i])*demPromedioTotal*pondNSC[0]/promInterpolN)
  interpolHoraS[i]=round((interpolHoraS[i])*demPromedioTotal*pondNSC[1]/promInterpolS)
  interpolHoraC[i]=round((interpolHoraC[i])*demPromedioTotal*pondNSC[2]/promInterpolC)

listaD=interpol(interpolHoraN,interpolHoraS,interpolHoraC,n,var)

#Algoritmo
for i in range(len(listaD)):
  resultado.append(permutarALG(listaD[i],Q,p,C,r))

#Arreglos de informacion
objFunc, prodN, prodS, prodC, fluxNS, fluxSN, fluxSC, fluxCS, fluxNC, fluxCN, demN, demS, demC, demTotal=[], [], [], [], [], [], [], [], [], [], [], [], [], []

for i in range(len(listaD)):
  objFunc.append(resultado[i][0])
  prodN.append(resultado[i][1])
  prodS.append(resultado[i][2])
  prodC.append(resultado[i][3])
  fluxNS.append(resultado[i][4])
  fluxSN.append(resultado[i][5])
  fluxSC.append(resultado[i][6])
  fluxCS.append(resultado[i][7])
  fluxNC.append(resultado[i][8])
  fluxCN.append(resultado[i][9])
  demN.append(listaD[i][0])
  demS.append(listaD[i][1])
  demC.append(listaD[i][2])
  demTotal.append(sum(listaD[i]))

#Inicio grafos
X = np.arange(0, 24, 24/len(listaD))

#Datos importantes
fluxSum=[]
fluxSum.append(sum(fluxNS)*24/len(fluxNS))
fluxSum.append(sum(fluxSN)*24/len(fluxSN))
fluxSum.append(sum(fluxSC)*24/len(fluxSC))
fluxSum.append(sum(fluxCS)*24/len(fluxCS))
fluxSum.append(sum(fluxNC)*24/len(fluxNC))
fluxSum.append(sum(fluxCN)*24/len(fluxCN))

print("Costo total en el d\'ia: $",round(sum(objFunc)*24/len(listaD),2))
print("Flujo total:",round(sum(fluxSum),2),"[MW]")
print("Energia perdida en flujo:", round(sum(fluxSum)*r,2),"[MW]")
print("=====Desglose flujos=====")
print("Flujo total NS:",round(fluxSum[0],2),"[MW]")
print("Flujo total SN:",round(fluxSum[1],2),"[MW]")
print("Flujo total SC:",round(fluxSum[2],2),"[MW]")
print("Flujo total CS:",round(fluxSum[3],2),"[MW]")
print("Flujo total NC:",round(fluxSum[4],2),"[MW]")
print("Flujo total CN:",round(fluxSum[5],2),"[MW]")


#Grafo demandas
plt.plot(X, demN, color='r', label='demN')
plt.plot(X, demS, color='g', label='demS')
plt.plot(X, demC, color='b', label='demC')
plt.plot([0,24], [0,0], color='0', linestyle='dotted') 
plt.xlabel("Tiempo [h]") 
plt.ylabel("Demanda [MW/h]") 
plt.title("Demanda a lo largo del tiempo") 
plt.legend() 
plt.show()

#Grafo produccion en cada nodo
plt.plot(X, prodN, color='r', label='qN') 
plt.plot(X, prodS, color='g', label='qS')
plt.plot(X, prodC, color='b', label='qC')
plt.plot([0,24],[Q[0],Q[0]], linestyle='dotted', color='r')
plt.plot([0,24],[Q[1],Q[1]], linestyle='dotted', color='g')
plt.plot([0,24],[Q[2],Q[2]], linestyle='dotted', color='b')
plt.plot([0,24], [0,0], color='0', linestyle='dotted')
plt.xlabel("Tiempo [h]") 
plt.ylabel("Produccion [MW/h]") 
plt.title("Produccion en cada nodo") 
plt.legend() 
plt.show()

#Grafo Norte
plt.plot(X, demN, color='r', label='demN')
plt.plot(X, prodN, color='g', label='qN') 
plt.plot([0,24], [0,0], color='0', linestyle='dotted') 
plt.plot([0,24], [Q[0]*p[0][0],Q[0]*p[0][0]], color='c', linestyle='dotted', label='tSolar')
plt.plot([0,24], [Q[0]*(p[0][0]+p[0][1]),Q[0]*(p[0][0]+p[0][1])], color='m', linestyle='dotted', label='tCarbon') 
plt.plot([0,24], [Q[0],Q[0]], color='y', linestyle='dotted', label='tGas') 
plt.xlabel("Tiempo [h]") 
plt.ylabel("Energia [MW/h]") 
plt.title("Produccion y demanda en el Norte") 
plt.legend() 
plt.show()

#Grafo Sur
plt.plot(X, demS, color='r', label='demS')
plt.plot(X, prodS, color='g', label='qS') 
plt.plot([0,24], [0,0], color='0', linestyle='dotted')
plt.plot([0,24], [Q[1]*p[1][0],Q[1]*p[1][0]], color='c', linestyle='dotted', label='tSolar')
plt.plot([0,24], [Q[1]*(p[1][0]+p[1][1]),Q[1]*(p[1][0]+p[1][1])], color='m', linestyle='dotted', label='tCarbon') 
plt.plot([0,24], [Q[1],Q[1]], color='y', linestyle='dotted', label='tGas') 
plt.xlabel("Tiempo [h]") 
plt.ylabel("Energ\'ia [MW/h]") 
plt.title("Producc\'ion y demanda en el Sur") 
plt.legend() 
plt.show()

#Grafo Centro
plt.plot(X, demC, color='r', label='demC')
plt.plot(X, prodC, color='g', label='qC') 
plt.plot([0,24], [0,0], color='0', linestyle='dotted') 
plt.plot([0,24], [Q[2]*p[2][0],Q[2]*p[2][0]], color='c', linestyle='dotted', label='tSolar')
plt.plot([0,24], [Q[2]*(p[2][0]+p[2][1]),Q[2]*(p[2][0]+p[2][1])], color='m', linestyle='dotted', label='tCarbon') 
plt.plot([0,24], [Q[2],Q[2]], color='y', linestyle='dotted', label='tGas') 
plt.xlabel("Tiempo [h]") 
plt.ylabel("Energ\'ia [MW/h]") 
plt.title("Producc\'ion y demanda en el Centro") 
plt.legend() 
plt.show()

#Grafo flujos
plt.plot(X, fluxNS, color='r', label='NS')
plt.plot(X, fluxSC, color='b', label='SC')
plt.plot(X, fluxNC, color='g', label='NC')
plt.plot(X, fluxSN, color='c', label='SN')
plt.plot(X, fluxCS, color='m', label='CS')
plt.plot(X, fluxCN, color='y', label='CN') 
plt.plot([0,24], [0,0], color='0', linestyle='dotted') 
plt.xlabel("Tiempo [h]") 
plt.ylabel("Flujo [MW/h]") 
plt.title("Flujos a lo largo del tiempo") 
plt.legend() 
plt.show()

#Grafo funcion objetivo
plt.plot(X, objFunc)
plt.xlabel("Tiempo [h]") 
plt.ylabel("Costo instantaneo [$/h]") 
plt.title("Funcion objetivo")
plt.show()
