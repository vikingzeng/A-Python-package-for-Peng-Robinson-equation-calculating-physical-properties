#Copyright (c) 2025 Lau JH. All rights reserved.
#This script is used to calculate fugacity coefficients for the pure or the mixture.
#It's prohibitively to put it into any form of commercial use.\n#It's also highly not recommend to share with people other than Xufeng Lin's group members and alumni or people who has acquired Lau JH's permission

import os
import math
import numpy as np
import sympy as sp

R=8.31451

def is_float(num)->bool:
    try:
        float(num)
        return True
    except ValueError:
        return False
def is_int(num)->bool:
    try:
        int(num)
        return True
    except ValueError:
        return False
    
def kappa_calculator(omega)->float:
    return 0.37464+1.54226*omega-0.26992*omega**2
def alpha_calculator(kappa,Tr)->float:
    return pow((1+kappa*(1-math.sqrt(Tr))),2)
def a_Tc_calculator(R,Tc,Pc)->float:
    return 0.45724*pow(R,2)*pow(Tc,2)/Pc
def b_Tc_calculator(R,Tc,Pc)->float:
    return 0.07780*R*Tc/Pc
def a_T_calculator(a_Tc,alpha)->float:
    return a_Tc*alpha
def A_calculator(a_T,R,P,T)->float:
    return (a_T*P)/pow(R*T,2)
def B_calculator(b_T,R,P,T)->float:
    return (b_T*P)/(R*T)
def mi_calculator(omega)->float:
    return 0.480+1.574*omega-0.0176*pow(omega,2)

def deltaij_calculator(list_deltat,list_deltaj,list_deltaij):
    for i in range(0,n_component):
        for j in range(0,n_component):
            list_deltaij[i][j]=list_deltat[i][j]+list_deltaj[i][j]
        
def deltat_calculator(list_deltat,case):
    match case:
        case 1:
            for i in range(0,n_component):
                for j in range(0,n_component):
                    if j!=i:
                        list_deltat[i][j]=1-(list_a_case1[i]-list_b_case1[i]*(T/list_Tc[i]-1)+list_Gij[i][j])/(1+list_mi[i]*(1-pow(T/list_Tc[i],0.5)))
                    else:
                        list_deltat[i][j]=0
        case 2:
            for i in range(0,n_component):
                for j in range(0,n_component):
                    if j!=i:
                        list_deltat[i][j]=list_a1_case2[i]*pow(Tr,2)+list_b1_case2[i]*Tr+list_c1_case2[i]
                    else:
                        list_deltat[i][j]=0
                    
        case 3:
            for i in range(0,n_component):
                for j in range(0,n_component):
                    if j!=i:
                        list_deltat[i][j]=list_a_case3[i]*T+list_b_case3[i]
                    else:
                        list_deltat[i][j]=0
                    
            
            
def deltaj_calculator(list_deltaj,case,n_component):
    match case:
        case 1:
            return 
        case 2:
            for i in range(0,n_component):
                for j in range(0,n_component):
                    if j!=i:
                        list_deltaj[i][j]=list_a_case2[i]*pow(list_omega[i],2)+list_b_case2[i]*list_omega[i]+list_c_case2[i]
                    else:
                        list_deltat[i][j]=0
        case 3:
            return 
    
print("This script calculate fugacity coefficient according to Peng-Robinson equation of state")
print('You can cite: Peng, Ding-Yu, and Donald B. Robinson. "A new two-constant equation of state." Industrial & Engineering Chemistry Fundamentals 15.1 (1976): 59-64. If you need it.')
print('The universal gas constant is chosen to be 8.31451 J·mol-1·K-1')
print('The supporting type of material are:\n 1. pure \n 2. mixture')

list_Tc=[]
list_PC=[]
list_omega=[]
list_fraction=[]
list_kappa=[]
list_alpha=[]
list_a_Tc=[]
list_a_T=[]
list_b_Tc=[]
list_Ai=[]
list_Bi=[]
list_deltat=[]
list_deltaj=[]
list_deltaij=[]
list_a_case1=[]
list_b_case1=[]
list_mi=[]
list_Gij=[]
list_a_case2=[]
list_b_case2=[]
list_c_case2=[]
list_a1_case2=[]
list_b1_case2=[]
list_c1_case2=[]
list_a_case3=[]
list_b_case3=[]


T=0
P=0
Tr=0
n_component=1
case=0
sum_b=0
sum_a=0

def pure(list_Tc,list_PC,list_omega,list_fraction,n_component,sum_a,sum_b):
    flag_pure=n_component
    count_species=1
    flag_once=True
    while True:
        print('please type in the temperature in the unit of K')
        T=input()
        if is_float(T):
            T=float(T)
            break
        else:
            print(T+" is not a digit, please try again")
            continue  
    while True:
        print('please type in the pressure in the unit of Pa')
        P=input()
        if is_float(P):
            P=float(P)
            break
        else:
            print(P+" is not a digit, please try again")
            continue 
        
    while n_component>0:
        case=0
        while True and flag_once==True and n_component>1:
            print("Binary Interaction Parameters are computed according to 'Moysan, J. M.; Paradowski, H.; Vidal, J. Prediction of phase behaviour of gas-containing systems with cubic equations of state. Chem. Eng. Sci. 1986, 41, 2069−2074'")
            print("please look up parameters you need in this article.")
            print('please choose which case your system is:\n 1. System with H2, CH4, N2 and CO'+
                  '\n 2. System with CO2 and H2S\n 3. System with H2O or methanol \n\notherwise, for systems with weak solvent effect type 1 is suggested')
            case=input()
            if is_int(case) and int(case) in (1,2,3):
                flag_once=False
                break
            else:
                print('your choice is not supported, please try again')
                continue
        while True:
            print('please type in the critical temperature of your '+str(count_species)+'th component in the unit of K')
            Tc=input()
            if  is_float(Tc):
                Tc=float(Tc)
                list_Tc.append(Tc)
                break
            else:
                print(Tc+" is not a digit, please try again")
                continue
        while True:
            print('please type in the critical pressure of your '+str(count_species)+'th component in the unit of Pa')
            Pc=input()
            if is_float(Pc):
                Pc=float(Pc)
                list_PC.append(Pc)
                break
            else:
                print(Pc+" is not a digit, please try again")
                continue
        while True:
            print('please type in the acentric factor of your '+str(count_species)+'th component')
            omega=input()
            if is_float(omega):
                omega=float(omega)
                list_omega.append(omega)
                break
            else:
                print(omega+" is not a digit, please try again")
                continue
        while True and flag_pure>1:
            print('please type in the mole ratio of your '+str(count_species)+'th component:')
            fraction=input()
            if is_float(fraction) and float(fraction)<1:
                fraction=float(fraction)
                list_fraction.append(fraction)
                break
            else:
                print(fraction+" is not a float or it is larger than 1, please try again")
                continue
        n_component=n_component-1
        count_species+=1
    
    n_component=flag_pure
    
    if case==1 and n_component>1:
        for i in range(0,n_component):
            print('please type in coefficient a for your '+str(i)+'th component')
            while True:
                a=input()
                if is_float(a):
                    a=float(a)
                    list_a_case1[i]=a 
                    break
                else:
                    print(a+'is not a float, please try again')
            print('please type in coefficient b for your '+str(i)+'th component')
            while True:
                b=input()
                if is_float(b):
                    b=float(b)
                    list_b_case1[i]=b
                    break
                else:
                    print(b+'is not a float, please try again')
        for i in range(0,n_component):
            for j in range(0,n_component) and j!=i:
                print('please type in the corrective term for your '+str(i)+'th component with '+str(j)+'th component')
                while True:
                    Gij=input()
                    if is_float(Gij):
                        Gij=float(Gij)
                        list_Gij[i].append(Gij)
                        break
                    else:
                        print(Gij+'is not a float, please try again')
                        continue
    if case==2 and n_component>1:
        for i in range(0,n_component):
            print('please type in coefficient a for the temperature effect on your '+str(i)+'th component')
            while True:
                a=input()
                if is_float(a):
                    a=float(a)
                    list_a_case2[i]=a 
                    break
                else:
                    print(a+'is not a float, please try again')
            print('please type in coefficient b for the temperature effect on your '+str(i)+'th component')
            while True:
                b=input()
                if is_float(b):
                    b=float(b)
                    list_b_case2[i]=b
                    break
                else:
                    print(b+'is not a float, please try again')
            print('please type in coefficient c for the temperature effect on your '+str(i)+'th component')
            while True:
                c=input()
                if is_float(c):
                    c=float(c)
                    list_c_case2[i]=c 
                    break
                else:
                    print(c+'is not a float, please try again')
            print('please type in coefficient a1 for the solvent effect on your '+str(i)+'th component')
            while True:
                a1=input()
                if is_float(a1):
                    a1=float(a1)
                    list_a1_case2[i]=a1 
                    break
                else:
                    print(a1+'is not a float, please try again')
            print('please type in coefficient b1 for the solvent effect on your '+str(i)+'th component')
            while True:
                b1=input()
                if is_float(b1):
                    b1=float(b1)
                    list_b1_case2[i]=b1 
                    break
                else:
                    print(b1+'is not a float, please try again')
            print('please type in coefficient c1 for the solvent effect on your '+str(i)+'th component')
            while True:
                c1=input()
                if is_float(c1):
                    c1=float(c1)
                    list_c1_case2[i]=c1 
                    break
                else:
                    print(c1+'is not a float, please try again')
    if case==3 and n_component>1:
        for i in range(0,n_component):
            print('please type in coefficient a2 for your '+str(i)+'th component')
            while True:
                a2=input()
                if is_float(a2):
                    a2=float(a2)
                    list_a_case3[i]=a2 
                    break
                else:
                    print(a2+'is not a float, please try again')
            print('please type in coefficient b2 for your '+str(i)+'th component')
            while True:
                b2=input()
                if is_float(b2):
                    b2=float(b2)
                    list_b_case3[i]=b2
                    break
                else:
                    print(b2+'is not a float, please try again')

    Tr=T/Tc
    if flag_pure==1:
        kappa=kappa_calculator(omega)
        alpha=alpha_calculator(kappa,Tr)
        a_Tc=a_Tc_calculator(R,Tc,Pc)
        b_Tc=b_Tc_calculator(R,Tc,Pc)
        a_T=a_T_calculator(a_Tc,alpha)
        b_T=b_Tc
        A=A_calculator(a_T,R,P,T)
        B=B_calculator(b_T,R,P,T)
        coeff=[
            1,
            -(1-B),
            (A-3*B**2-2*B),
            -(A*B-B**2-B**3)
            ]
        roots=np.roots(coeff)
        real_roots=[root.real for root in roots if abs(root.imag)<1e-10]
        if not real_roots:
            raise ValueError("There is no real root, please check your input")
        z=max(real_roots)
        if z - B > 0 and (z + (1 + math.sqrt(2) ) * B) / (z - (1 - math.sqrt(2) )* B) > 0: 
            ln_fi = (z - 1) - np.log(z  - B) - (A / (math.sqrt(2)  * 2 * B)) * np.log((z  + (1 + math.sqrt(2) )* B) / (z - (1 - math.sqrt(2)  )* B)) 
            fi = np.exp(ln_fi)  
            print(f"Fugacity coefficient: {fi}") 
        else: 
            print("Invalid input leads to negative argument in logarithm function. Please check your input.") 
        os.system('pause')

    if flag_pure>1:
        for i in range(0,n_component):
            if case==1:
                list_mi.append(mi_calculator(list_omega[i]))
            list_kappa.append(kappa_calculator(list_omega[i]))
            list_alpha.append(alpha_calculator(list_kappa[i],Tr))
            list_a_Tc.append(a_Tc_calculator(R,list_Tc[i],list_PC[i]))
            list_b_Tc.append(b_Tc_calculator(R,list_Tc[i],list_PC[i]))
            list_a_T.append(a_T_calculator(list_a_Tc[i],list_alpha[i]))
            list_Ai.append(A_calculator(list_a_T[i],R,P,T))
            list_Bi.append(B_calculator(list_b_Tc[i],R,P,T))
        for i in range(0,n_component):
            sum_b+=list_fraction[i]*list_b_Tc[i]
            for j in range(0,n_component):
                sum_a+=list_fraction[i]*list_fraction[j]*(1-list_deltaij[i][j])*pow(list_a_Tc[i],0.5)*pow(list_a_Tc[j],0.5)
        for i in range(0,n_component):
            A=list_Ai[i]
            B=list_Bi[i]
            coeff=[
                1,
                -(1-B),
                (A-3*B**2-2*B),
                -(A*B-B**2-B**3)
            ]
            roots=np.roots(coeff)
            real_roots=[root.real for root in roots if abs(root.imag)<1e-10]
            if not real_roots:
                raise ValueError("There is no real root, please check your input")
            z=max(real_roots)
            if z - B > 0 and (z + (1 + math.sqrt(2) ) * B) / (z - (1 - math.sqrt(2) )* B) > 0:
                sum_temp=0
                for j in range(n_component):
                    sum_temp=sum_temp+list_fraction[j]*(1-list_deltaij[j][i])*pow(list_a_Tc[j],0.5)*pow(list_a_Tc[i],0.5)
                ln_fi = (list_b_Tc[i]/sum_b) * (z - 1) - np.log(z  - B) - (A / (math.sqrt(2)  * 2 * B)) * (2*(sum_temp/sum_a)-(list_b_Tc[i]/sum_b)) *np.log((z  + (1 + math.sqrt(2) )* B) / (z - (1 - math.sqrt(2)  )* B)) 
                fi = np.exp(ln_fi)  
                print(f"Fugacity coefficient of {i}th component is : {fi} \n")   
            else: 
                print("Invalid input leads to negative argument in logarithm function. Please check your input.")  
        os.system('pause') 

choice=''
while True:
    print('please type in your choice:')
    choice=input()
    if int(choice)==1:
        pure(list_Tc,list_PC,list_omega,list_fraction,n_component,sum_a,sum_b)
        break
    elif int(choice)==2:
        print('please type in the number of your components:')
        flag=True
        while flag:
            temp=input()
            if is_int(temp):
                n_component=int(temp)
                for i in range(0,n_component):
                    list_temp=[0]*n_component
                    list_deltaij.append(list_temp)
                    list_deltat.append(list_temp)
                    list_deltaj.append(list_temp)
                    list_Gij.append(list_temp)
                pure(list_Tc,list_PC,list_omega,list_fraction,n_component,sum_a,sum_b)
                flag=False
            else:
                print(temp+" is not a digit, please try again")
                continue
        break
    else:
        print('your choice is not supported, please try again')