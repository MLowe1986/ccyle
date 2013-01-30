# import the pysb module and all its methods and functions
from pysb import *
from pysb.macros import catalyze, degrade
    
# instantiate a model
Model()


#Declare Monomers
Monomer('GF',['b'])                               #Growth Factor that activates synthesis of AP1
Monomer('AP1', ['b', 'y'],{'y':('i','a')} )      #Transcription factor that promotes the synthesis of cyclin D.
Monomer('Cd', ['b', 'y'], {'y':('i', 'a')})       #Cyclin D

#Declare Parameters
Parameter('kfap1', 1.663e-05)	    #Calculated from the Michaelis constant, Kagf of 0.1microMoles provided in the paper (Km = (kr + kcat)/kf)
Parameter('krap1', 1.0e-03)       #Assumed typical reverse rate constant
Parameter('kcap1', 1.0)           #Assumed typical catalytic rate constant 
#Parameter('kdap1', 4.17e-05)      #Apparent first-order rate constant for AP1 transcription factor degradation in h^-1

# initial conditions
Parameter('GF_0', 1000)      #the initial concentration of GF in number of molecules/pL - Converted from 1microM in paper to be the number of molecules in a cell (which is about 1 picoLiter)
Initial(GF(b=None), GF_0)
Parameter('AP1_0',100000)      #the initial concentration of the substrate, inactive AP1, in molecules/pL. Converted from the rate of synthesis of the transcription factor AP1 given in microM per hour and depending on growth factor GF: V(t) = (kcap1xGF_0xAP1_0) divided by (AP1_0+Km)	                              
Initial(AP1(b=None, y='i'),AP1_0)

Observable('obsAP1i', AP1(b=None, y='i'))    #observe the inactive form of AP1
Observable('obsAP1a', AP1(b=None, y='a'))    #observe the active form of AP1
Observable('obsGF', GF(b=None))             #observe Growth factor

#Catalyze Rule for reaction: The GF activates the synthesis of the transcription factor AP1. 
catalyze(GF(), 'b', AP1(y='i'), 'b', AP1(y='a'), (kfap1, krap1, kcap1))  # calls the function catalyze to simulate synthesis of AP1

#Rule for degradation of AP1 
#degrade(AP1(y='i'), kdap1)