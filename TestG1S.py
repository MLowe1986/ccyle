
from pysb import *
from pysb.macros import catalyze, degrade, equilibrate, bind, bind_table, synthesize
    
# instantiate a model
Model()

#Declare Monomers
Monomer('pRB', ['b', 'S'],{'S':['u', 'p', 'pp']})  # Retinoblastoma protein unphosphorylated, phosphorylated once, and phosphorylated twice. 
Monomer('Ce',['b', 'S'], {'S':['i', 'a']})
Monomer('Cdk2',['b', 'c', 'S'], {'S':['i', 'a']})
Monomer('p27_p21',['b', 'S'], {'S':['u', 'p']})
Monomer('Wee1',['b', 'S'], {'S':['a', 'i']})
Monomer('CDC25', ['b', 'S'], {'S':['a', 'i']})

#Declare Parameters

#Parameter for Synthesis of Cyclin E 
Parameter('kSCe', 7.5e-02)	    #This is k5 in Iwamoto paper

#Parameter for Degradation of Cyclin E
Parameter('kDCe', 2.5e-03)          #This is k6 in Iwamoto paper

#Cyclin E/Cdk2 binding parameters: 
Parameter('kcdk2f',1.25e-03)
Parameter('kcdk2r', 2.5e-04)  #forward and reverse rates respectively for formation of CyclinE/Cdk2 complex (k7 and k8 from Iwamoto- CHECK THIS BY COMPARING TO EQUATIONS IN PAPER!!)

#CyclinE/Cdk2 Complex Equilibrate Parameters:
###I AM NOT SURE WHAT THESE ARE IF THEY ARE CORRECT
Parameter('kfcecdk', 2.5e-02)     #k22 from Iwamoto = NOT SURE
Parameter('krcecdk', 1.75e-03)    #k23 from Iwamoto = NOT SURE

#Cyclin E/Cdk2/p27_p21 binding parameters:
#p27_rates = [2.5e-02, 1.75e-03]  #forward and reverse rates respectively for formation of CyclinE/Cdk2/p27_p21 complex (k22 and k23 from Iwamoto- CHECK THIS BY COMPARING TO EQUATIONS IN PAPER!!)

Parameter('kfcp27_p21', 1.0e-02)     #k35 from Iwamoto
Parameter('krcp27_p21', 1.74e-04)    #k25 from Iwamoto

#The Iwamoto paper seems to show a synthesis of p27_p21 instead of a phosphorylation. Which do I show? Or do I show the synthesis and degradation rates in Iwamoto as the forward/reverse rates for phosphorylation... this is what I have done below: 
Parameter('kfp27_p21', 5.0e-08)           #k34 from Iwamoto- NOT SURE
Parameter('krp27_p21', 2.5e-03 )          #k30 from Iwamoto - NOT SURE

#Paramerers for CDC25 equilibration:
Parameter('kfCDC25', 5.0e-08)                   #THESE PARAMETERS ARE NOT CORRECT- ESTIMATED FOR NOW. CDC25 IS NOT IN IWAMOTO PAPER
Parameter('krCDC25', 2.5e-03)                   #THESE PARAMETERS ARE NOT CORRECT- ESTIMATED FOR NOW. CDC25 IS NOT IN IWAMOTO PAPER



#Observables
#Observable('obsCdh1i', Cdh1(b=None, S='i'))    
#Observable('obsCdh1a', Cdh1(b=None, S='a'))    
#Observable('obsSkp2', Skp2(b=None))            
Observable('obsCei', Ce(b=None, S='i'))
Observable('obsCea', Ce(b=None, S='a'))
#Observable('obsE2Fu', E2F(b=None, S='u'))
#Observable('obsE2Fp', E2F(b=None, S='p'))
Observable('obspRBu', pRB(b=None, S='u'))
Observable('obspRBp', pRB(b=None, S='p'))
Observable('obspRBpp', pRB(b=None, S='pp'))
Observable('obsCdk2i', Cdk2(b=None, S='i'))
Observable('obsCdk2a', Cdk2(b=None, S='a'))
Observable('obsCe_Cdk2i', Ce(b=1, S='a') % Cdk2(b=1, S ='i'))
Observable('obsCe_Cdk2a', Ce(b=1, S='a') % Cdk2(b=1, S ='a'))
Observable('obsCe_Cdk2a_p27', Ce(b=1, S='a') % Cdk2(b=1, c=2, S ='a') % p27_p21(b=2, S='u'))
Observable('obsp27_p21u', p27_p21(b=None, S='u'))
Observable('obsp27_p21p', p27_p21(b=None, S='p'))
Observable('obsWee1i', Wee1(b=None, S='i'))
Observable('obsWee1a', Wee1(b=None, S='a'))
Observable('obsCDC25i', CDC25(b=None, S='i'))    
Observable('obsCDC25a', CDC25(b=None, S='a'))

#Initial Conditions

Parameter('Ce_0', 1.0e-03)           #Y1 in Iwamoto
Initial(Ce(b=None, S='a'), Ce_0)

Parameter('Cdk2_0', 1.35e01)           #Y4 in Iwamoto
Initial(Cdk2(b=None, c=None, S='i'),Cdk2_0)   

Parameter('Ce_Cdk2_0', 1.0e-03)                       #Y6 in Iwamoto
Initial(Ce(b=1, S='a') % Cdk2(b=1, c=None, S='a'), Ce_Cdk2_0)   

Parameter('p27_21_0', 1.0e01)                   #Y10 in Iwamoto
Initial(p27_p21(b=None, S='u'),p27_21_0)

Parameter('Ce_Cdk2_p27_p21_0', 1.0e00)                   #Y12 in Iwamoto
Initial(Ce(b=1, S='a') % Cdk2(b=1, c=2, S='a') % p27_p21(b=2, S='u'), Ce_Cdk2_p27_p21_0)

Parameter('CDC25_0', 1.0e01)				#THIS PARAMETER IS NOT CORRECT- ESTIMATED FOR NOW- NOT IN IWAMOTO PAPER	
Initial(CDC25(b=None, S='i'), CDC25_0)     

#Declare Rules

#E2F promotes the synthesis of Cyclin E  
#I can show synthesis, as shown below, but how do I show the promotion part????

#Cyclin E Synthesis and Degradation
synthesize(Ce(b=None, S = 'a'), kSCe)
degrade(Ce(b=None, S='a'),kDCe)

#Cyclin E binds to Cdk2 
Rule('Ce_Cdk2_bind',Ce(b=None, S='a') + Cdk2(b=None, c=None, S='i') <> Ce(b=1, S='a') % Cdk2(b=1, c=None, S ='i'), *[kcdk2f, kcdk2r])

#The Ce/Cdk2 complex equilibrates from inactive to active
equilibrate(Ce(b=1, S='a') % Cdk2(b=1, c=None, S='i'),Ce(b=1, S='a') % Cdk2(b=1, c=None, S='a'),[kfcecdk,krcecdk])

#The Ce/Cdk2 complex binds to p27_p21
Rule ('Ce_Cdk2_p27_p21_bind',p27_p21(b=None, S='u') + Ce(b=1, S='a') % Cdk2(b=1, c=None, S='a') <> Ce(b=1, S='a') % Cdk2(b=1, c=2, S ='a') % p27_p21(b=2, S='u') , *[kfcp27_p21, krcp27_p21])

#Phosphorylate p27_21
equilibrate(p27_p21(b=None,S='u'),p27_p21(b=None, S='p'),[kfp27_p21, krp27_p21])

#Equilibrate CDC25
equilibrate(CDC25(b=None,S='i'),CDC25(b=None, S='a'),[kfCDC25, krCDC25])

#Previously, I thought synthesis should be shown as an enzyme reaction because it is promoted by a certain factor:
#catalyze(E2F(S='u'), 'b', Ce(S='i'), 'b', Ce(S='a'), (kfce, krce, kcce))  # calls the function catalyze to simulate activation of Cyclin E

#Parameters that bind to Ce
#bind_table([[ Ce(S= 'a')]
#	     [pRB(S = 'u'), pRB_rates, None]              #pRB inhibits Cyclin E 
#	     [pRB(S = 'p'), pRBp_rates, None]             #pRBp inhibits Cyclin E 
#	     [cdk2(S = 'i'), None, cdk2_rates]])	  #cdk2 binds to Cylcin E to form a complex. 