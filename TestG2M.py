from pysb import *
from pysb.macros import catalyze, degrade, equilibrate, bind, bind_table, synthesize
    
# instantiate a model
Model()

#Declare Monomers
Monomer('Cb',['b', 'S'], {'S':['i', 'a']})
Monomer('Cdk1',['b', 'c', 'S'], {'S':['i', 'a']})
Monomer('p27_p21',['b', 'S'], {'S':['u', 'p']})
Monomer('Wee1',['b', 'S'], {'S':['a', 'i']})
Monomer('CDC25', ['b', 'S'], {'S':['a', 'i']})
Monomer('CDC20', ['b', 'S'], {'S':['a', 'i']})
Monomer('Cdh1', ['b', 'S'], {'S':['a', 'i']})


#Declare Parameters

#Parameter for Synthesis of Cyclin B
Parameter('kSCb', 8.0e-04)	    #FILL IN VALUE- THIS VALUE IS NOT IN IWAMOTO PAPER

#Parameter for Degradation of Cyclin B
Parameter('kDCb', 5.0e-04)          #FILL IN VALUE- THIS VALUE IS NOT IN IWAMOTO PAPER

#Cyclin B/Cdk1 binding parameters: 
Parameter('kcdk1f',1.0e-03)
Parameter('kcdk1r', 2.0e-04)  #FILL IN VALUES- THESE VALUES ARE NOT IN IWAMOTO PAPER

#CyclinB/Cdk1 Complex Equilibrate Parameters:

Parameter('kfcbcdk', 1.9e-02)     #FILL IN VALUE- THIS VALUE IS NOT IN IWAMOTO PAPER
Parameter('krcbcdk', 5.0e-04)    #FILL IN VALUE- THIS VALUE IS NOT IN IWAMOTO PAPER

#Cyclin B/Cdk1/p27_p21 binding parameters:
Parameter('kfcp27_p21', 1.5e-03)     #FILL IN VALUE- THIS VALUE IS NOT IN IWAMOTO PAPER
Parameter('krcp27_p21', 1.75e-04)    #FILL IN VALUE- THIS VALUE IS NOT IN IWAMOTO PAPER

#Paramerers for CDC25 equilibration:
Parameter('kfCDC25', 5.0e-08)                   #THESE PARAMETERS ARE NOT CORRECT- ESTIMATED FOR NOW. CDC25 IS NOT IN IWAMOTO PAPER
Parameter('krCDC25', 2.5e-03)                   #THESE PARAMETERS ARE NOT CORRECT- ESTIMATED FOR NOW. CDC25 IS NOT IN IWAMOTO PAPER

#Parameters for CDC20 equilibration:
Parameter('kfCDC20', 5.0e-08)                   #THESE PARAMETERS ARE NOT CORRECT- ESTIMATED FOR NOW. CDC20 IS NOT IN IWAMOTO PAPER
Parameter('krCDC20', 2.5e-03)                   #THESE PARAMETERS ARE NOT CORRECT- ESTIMATED FOR NOW. CDC20 IS NOT IN IWAMOTO PAPER

#Parameters for Wee1 equilibration:
Parameter('kfWee1', 5.0e-08)                   #THESE PARAMETERS ARE NOT CORRECT- ESTIMATED FOR NOW. Check Iwamoto
Parameter('krWee1', 2.5e-03)                   #THESE PARAMETERS ARE NOT CORRECT- ESTIMATED FOR NOW. Check Iwamoto



#Observables
#Observable('obsCdh1i', Cdh1(b=None, S='i'))    
#Observable('obsCdh1a', Cdh1(b=None, S='a'))               
Observable('obsCbi', Cb(b=None, S='i'))
Observable('obsCba', Cb(b=None, S='a'))
Observable('obsCb_Cdk1i', Cb(b=1, S='a') % Cdk1(b=1, S ='i'))
Observable('obsCb_Cdk1a', Cb(b=1, S='a') % Cdk1(b=1, S ='a'))
Observable('obsCb_Cdk1a_p27', Cb(b=1, S='i') % Cdk1(b=1, c=2, S ='i') % p27_p21(b=2, S='u'))
Observable('obsp27_p21u', p27_p21(b=None, S='u'))
Observable('obsp27_p21p', p27_p21(b=None, S='p'))
Observable('obsWee1i', Wee1(b=None, S='i'))
Observable('obsWee1a', Wee1(b=None, S='a'))
Observable('obsCDC25i', CDC25(b=None, S='i'))    
Observable('obsCDC25a', CDC25(b=None, S='a'))
Observable('obsCDC20i', CDC25(b=None, S='i'))    
Observable('obsCDC20a', CDC25(b=None, S='a'))
Observable('obsCdk1i', Cdk1(b=None, c= None, S='i'))
Observable('obsCdk1a', Cdk1(b=None, c= None, S='a'))


#Initial Conditions

Parameter('Cb_0', 4.0e-05)           #VALUE NOT CORRECT
Initial(Cb(b=None, S='a'), Cb_0)

Parameter('Cdk1_0', 1.35e01)           #VALUE NOT CORRECT
Initial(Cdk1(b=None, c=None, S='i'),Cdk1_0)   

Parameter('Cb_Cdk1i_0', 4.0e-04)                       #VALUE NOT CORRECT
Initial(Cb(b=1, S='a') % Cdk1(b=1, c=None, S='i'), Cb_Cdk1i_0)   

Parameter('Cb_Cdk1a_0', 1.0e-04)                       ##VALUE NOT CORRECT
Initial(Cb(b=1, S='a') % Cdk1(b=1, c=None, S='a'), Cb_Cdk1a_0)   

Parameter('p27_21_0', 1.0e01)                   #Y10 in Iwamoto
Initial(p27_p21(b=None, S='u'),p27_21_0)

Parameter('Cb_Cdk1_p27_p21_0', 1.0e-04)                   #VALUE NOT CORRECT
Initial(Cb(b=1, S='a') % Cdk1(b=1, c=2, S='a') % p27_p21(b=2, S='u'), Cb_Cdk1_p27_p21_0)

Parameter('CDC25_0', 1.0e01)				#THIS PARAMETER IS NOT CORRECT- ESTIMATED FOR NOW- NOT IN IWAMOTO PAPER	
Initial(CDC25(b=None, S='i'), CDC25_0)     

#Declare Rules


#Cyclin B Synthesis and Degradation
synthesize(Cb(b=None, S = 'a'), kSCb)
degrade(Cb(b=None, S='a'),kDCb)

#Cyclin B binds to Cdk1
Rule('Cb_Cdk1_bind',Cb(b=None, S='a') + Cdk1(b=None, c=None, S='i') <> Cb(b=1, S='a') % Cdk1(b=1, c=None, S ='i'), *[kcdk1f, kcdk1r])

#The Cb/Cdk1 complex equilibrates from inactive to active
equilibrate(Cb(b=1, S='a') % Cdk1(b=1, c=None, S='i'),Cb(b=1, S='a') % Cdk1(b=1, c=None, S='a'),[kfcbcdk,krcbcdk])

#The Cb/Cdk1 complex binds to p27_p21 - this complex is inactive
Rule ('Cb_Cdk1_p27_p21_bind',p27_p21(b=None, S='u') + Cb(b=1, S='a') % Cdk1(b=1, c=None, S='a') <> Cb(b=1, S='i') % Cdk1(b=1, c=2, S ='i') % p27_p21(b=2, S='u') , *[kfcp27_p21, krcp27_p21])

#Equilibrate CDC25
equilibrate(CDC25(b=None,S='i'),CDC25(b=None, S='a'),[kfCDC25, krCDC25])

#Equilibrate CDC20
equilibrate(CDC20(b=None,S='i'),CDC20(b=None, S='a'),[kfCDC20, krCDC20])

#Equilibrate Wee1
equilibrate(Wee1(b=None,S='i'),Wee1(b=None, S='a'),[kfWee1, krWee1])

##JUST NEED TO SHOW THE EFFECTS OF MONOMERS AND COMPLEXES ON OTHER MONOMERS AND COMPLEXES NOW