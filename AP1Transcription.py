# import the pysb module and all its methods and functions
from pysb import *
from pysb.macros import catalyze, degrade
    
# instantiate a model
Model()

#Declare Monomers
Monomer('GF',['b'])        #Growth Factor that activates synthesis of AP1
Monomer('AP1', ['b', 'y'],{'y':('i','a')} )      #Transcription factor that promotes the synthesis of cyclin D.

#Declare Parameters
Parameter('kf', 1.0e07)	    #Calculated from the Michaelis constant, Kagf of 0.1microMoles provided in the paper (Km = (kr + kcat)/kf)
Parameter('kr', 1.0e-03)    #Assumed typical reverse rate constant
Parameter('kc', 1.0)	    #Assumed typical catalytic rate constant 
Parameter('kdap1', 0.15)    #Apparent first-order rate constant for AP1 transcription factor degradation in h^-1


# initial conditions
Parameter('GF_0', 602200)      #the initial concentration of GF is 1microMole
Initial(GF(b=None), GF_0)
Parameter('AP1_0', 10000)           #WE NEED THIS VALUE!!!
Initial(AP1(b=None, y='i'),AP1_0)

#Catalyze Rule for reaction: The GF activates the synthesis of the transcription factor AP1. 
catalyze(GF(), 'b', AP1(y='i'), 'b', AP1(y='a'), (kf, kr, kc))  # calls the function catalyze to simulate synthesis of AP1
  
#Rule for degradation of AP1 
#degrade(AP1(), kdap1)
