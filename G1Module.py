# import the pysb module and all its methods and functions
from pysb import *
from pysb.macros import catalyze, degrade
    
# instantiate a model
Model()


#Declare Monomers
Monomer('GF',['b'])                               #Growth Factor that activates synthesis of AP1
Monomer('AP1', ['b', 'S'],{'S':['i','a']})       #Transcription factor that promotes the synthesis of cyclin D.
Monomer('Cd', ['b', 'S'], {'S':['i', 'a']})       #Cyclin D
Monomer('E2F', ['b', 'S'], {'S': ['u', 'p']})     #transcription factor E2F and the phosphorylated form of E2F
Monomer('pRB', ['b', 'S'],{'S':['u', 'p', 'pp']})  # Retinoblastoma protein unphosphorylated, phosphorylated once, and phosphorylated twice. 
Monomer('Cdk4_6', ['b'], {'S':['a', 'i']})          #Cdk4_6 are treated as one monomer



#AP1 Activation Parameters 
Parameter('kfap1', 1.663e-05)	    #Calculated from the Michaelis constant, Kagf of 0.1microMoles provided in the paper (Km = (kr + kcat)/kf)
Parameter('krap1', 1.0e-03)       #Assumed typical reverse rate constant
Parameter('kcap1', 1.0)           #Assumed typical catalytic rate constant 
Parameter('kdap1', 4.17e-05)      #Apparent first-order rate constant for AP1 transcription factor degradation in s^-1

#CyclinD Activation by AP1 Parameters
Parameter('kfcd1', 1.0e07)	    #Assumed 
Parameter('krcd1', 1.0e-03)    #Assumed typical reverse rate constant
Parameter('kccd1', 0.4)	    #Rate constant for Cyclin D synthesis induced by AP1 in h^-1  
Parameter('kddd1', 0.005)    #Apparent first-order rate constant for non-specific cyclin D protein degradation in s^-1

#Cyclin D Activation Parameters by E2F
Parameter('kfcd2', 1.0e07)	    #Assumed 
Parameter('krcd2', 1.0e-03)    #Assumed typical reverse rate constant
Parameter('kccd2', 0.4)	    #Rate constant for Cyclin D synthesis induced by AP1 in h^-1  
Parameter('kddd2', 0.005)    #Apparent first-order rate constant for non-specific cyclin D protein degradation in s^-1

#Cyclin D Inhibition Parameters
Parameter('kfi7', 60.2)   #forward rate for CyclinD inhibition by pRB    
Parameter('kri7', 1.0e-03)    #reverse rate for CyclinD inhibition by pRB
Parameter('kfi8', 1204.4)	    #forward rate for CyclinD inhibition by pRBp
Parameter('kri8', 1.0e-03)    #reverse rate for CyclinD inhibition by pRBp

#CyclinD/Cdk4_6 Parameters
Parameter('kcom1', )         #forward rate for CyclinD/Cdk4_6 complex formation   
Parameter('kdecom1', )    #reverse rate for CyclinD/Cdk4_6 complex formation   

# AP1 initial conditions
Parameter('GF_0', 602200)      #the initial concentration of GF in number of molecules/pL - Converted from 1microM in paper to be the number of molecules in a cell (which is about 1 picoLiter)
Initial(GF(b=None), GF_0)
Parameter('AP1_0',16.73)      #the initial concentration of the substrate, inactive AP1, in molecules/pL. Converted from the rate of synthesis of the transcription factor AP1 given in microM per hour and depending on growth factor GF: V(t) = (kcap1xGF_0xAP1_0) divided by (AP1_0+Km)	                              
Initial(AP1(b=None, S='i'),AP1_0)

#Observables
Observable('obsAP1i', AP1(b=None, S='i'))    #observe the inactive form of AP1
Observable('obsAP1a', AP1(b=None, S='a'))    #observe the active form of AP1
Observable('obsGF', GF(b=None))             #observe Growth factor
Observable('obsCdi', Cd(b=None, S='i'))
Observable('obsCda', Cd(b=None, S='a'))
Observable('obsE2Fu', E2F(b=None, S='u'))
Observable('obsE2Fp', E2F(b=None, S='p'))
Observable('obspRBu', pRB(b=None, S='u'))
Observable('obspRBp', pRB(b=None, S='p'))
Observable('obspRBpp', pRB(b=None, S='pp'))



#Catalyze Rule for reaction: The GF activates the synthesis of the transcription factor AP1. 
catalyze(GF(), 'b', AP1(S='i'), 'b', AP1(S='a'), (kfap1, krap1, kcap1))  # calls the function catalyze to simulate synthesis of AP1

#AP1 activates Cyclin D.
catalyze(AP1(S='a'), 'b', Cd(S='i'), 'b', Cd(S='a'), (kfcd1, krcd1, kccd1))  # calls the function catalyze to simulate synthesis of AP1 
  
#E2F activates Cyclin D  
catalyze(E2F(S='u'), 'b', Cd(S='i'), 'b', Cd(S='a'), (kfcd2, krcd2, kccd2))  # calls the function catalyze to simulate synthesis of AP1   
  
#pRB inhibits Cyclin D 
Rule('pRB_Cd__bind', pRB(b=None, S= 'u') + Cd(b=None, S = 'i') <> pRB(b=1, S='u') % Cd(b=1, S= 'i'), *[kfi7, kri7])  
#pRBp inhibits Cyclin D 
Rule('pRBp_Cd__bind', pRB(b=None, S= 'p') + Cd(b=None, S = 'i') <> pRB(b=1, S='p') % Cd(b=1, S= 'i'), *[kfi8, kri8])  

Rule('Cd_Cdk4_6_bind', Cd(b=None) + Cdk4_6(b=None, S='i') <> Cd(b=1) % Cdk4_6(b=1, S= 'i'), *[kcom1, kdecom1])

Rule('Cd_Cdk4_6active', Cd(b=1) % Cdk4_6(b=1, S='i') <> Cd(b=1) % Cdk4_6(b=1, S ='a'), *[???????????])

Rule('Cd_Cdk4_6_p27_p21_bind', Cd(b=1) % Cdk4_6(b=1, S='a') <> Cd(b=1) % Cdk4_6(b=1, S ='a') % p27_p21(b= , *[???????????])
  
#Rule for degradation of AP1 
degrade(AP1(), kdap1)






