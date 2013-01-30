# import the pysb module and all its methods and functions
# import the pysb module and all its methods and functions
from pysb import *
from pysb import catalyze   

# instantiate a model
Model()

# declare monomers
#AP1 shouldn't have to be declared again because it was declared in previous program- otherwise it is inefficient. Maybe put together as sections with declared monomers first. Since they all affect eachother, we can't really do each one as a separate program anyway... 
Monomer('AP1', ['b', 'y'],{'y':('i','a')} )      #Transcription factor that promotes the synthesis of cyclin D.
Monomer('Cd', ['b', 'y'], {'y':('i', 'a')})       #Cyclin D



#Declare Parameters
#PARAMETERS WILL HAVE TO BE GIVEN DIFFERENT NAMES IF THESE MODULES ARE PUT TOGETHER
#STILL NEED TO COME UP WITH ASSUMED PARAMETER VALUES FOR kf AND kr. 
Parameter('kf', 1.0e07)	    #Assumed 
Parameter('kr', 1.0e-03)    #Assumed typical reverse rate constant
Parameter('kc', 0.4)	    #Rate constant for Cyclin D synthesis induced by AP1 in h^-1  
Parameter('kddd', 0.005)    #Apparent first-order rate constant for non-specific cyclin D protein degradation in h^-1


# initial conditions
#DO WE NEED TO DECLARE INITIAL CONDITIONS, SINCE AP1 SHOULD HAVE A VALUE??


#Catalyze Rule for reaction: The GF activates the synthesis of the transcription factor AP1. 
catalyze(AP1(y='a'), 'b', Cd(y='i'), 'b', Cd(y='a'), (kf, kr, kc))  # calls the function catalyze to simulate synthesis of AP1
  
#Rule for degradation of AP1 
degrade(Cd(), kddd)


