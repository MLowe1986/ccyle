from pysb import *
from pysb.macros import catalyze, degrade, equilibrate, bind, bind_table, synthesize
    
# instantiate a model
Model()

#Declare Monomers
Monomer('pRB', ['b', 'S'],{'S':['u', 'p', 'pp']})  # Retinoblastoma protein unphosphorylated, phosphorylated once, and phosphorylated twice. 
Monomer('Ca',['b', 'S'], {'S':['i', 'a']})
Monomer('Cdk2',['b', 'c', 'S'], {'S':['i', 'a']})
Monomer('p27_p21',['b', 'S'], {'S':['u', 'p']})
Monomer('Wee1',['b', 'S'], {'S':['a', 'i']})
Monomer('CDC25', ['b', 'S'], {'S':['a', 'i']})
Monomer('E2F', ['b'])
Monomer('CDC20', ['b', 'S'], {'S':['a', 'i']})
Monomer('Cdh1', ['b', 'S'], {'S':['a', 'i']})
Monomer('Cb',['b', 'S'], {'S':['i', 'a']})
Monomer('Cdk1',['b', 'c', 'S'], {'S':['i', 'a']})

#Declare Parameters

#Parameter for Synthesis of Cyclin A
Parameter('kSCa', 8.0e-04)	    #This is k9 in Iwamoto paper

#Parameter for Degradation of Cyclin A
Parameter('kDCa', 5.0e-04)          #This is k10 in Iwamoto paper

#Cyclin A/Cdk2 binding parameters: 
Parameter('kcdk2f',1.0e-03)
Parameter('kcdk2r', 2.0e-04)  #forward and reverse rates respectively for formation of CyclinA/Cdk2 complex (k11 and k12 from Iwamoto- CHECK THIS BY COMPARING TO EQUATIONS IN PAPER!!)

#CyclinA/Cdk2 Complex Equilibrate Parameters:
###I AM NOT SURE WHAT THESE ARE IF THEY ARE CORRECT
Parameter('kfcacdk', 1.9e-02)     #k28 from Iwamoto = NOT SURE
Parameter('krcacdk', 5.0e-04)    #k29 from Iwamoto = NOT SURE

#Cyclin A/Cdk2/p27_p21 binding parameters:
Parameter('kfcp27_p21', 1.5e-03)     #k36 from Iwamoto
Parameter('krcp27_p21', 1.75e-04)    #k31 from Iwamoto

#The Iwamoto paper seems to show a synthesis of p27_p21 instead of a phosphorylation. Which do I show? Or do I show the synthesis and degradation rates in Iwamoto as the forward/reverse rates for phosphorylation... this is what I have done below: 
Parameter('kfp27_p21', 5.0e-08)           #k34 from Iwamoto- NOT SURE
Parameter('krp27_p21', 2.5e-03 )          #k30 from Iwamoto - NOT SURE

#Paramerers for CDC25 equilibration:
Parameter('kfCDC25', 5.0e-08)                   #THESE PARAMETERS ARE NOT CORRECT- ESTIMATED FOR NOW. CDC25 IS NOT IN IWAMOTO PAPER
Parameter('krCDC25', 2.5e-03)                   #THESE PARAMETERS ARE NOT CORRECT- ESTIMATED FOR NOW. CDC25 IS NOT IN IWAMOTO PAPER

#Parameters for Cdh1 equilibration: 
Parameter('kfCdh1', 5.0e-08)                   #THESE PARAMETERS ARE NOT CORRECT- ESTIMATED FOR NOW. Cdh1 IS NOT IN IWAMOTO PAPER
Parameter('krCdh1', 2.5e-03)                   #THESE PARAMETERS ARE NOT CORRECT- ESTIMATED FOR NOW. Cdh1 IS NOT IN IWAMOTO PAPER

#Observables
#Observable('obsCdh1i', Cdh1(b=None, S='i'))    
#Observable('obsCdh1a', Cdh1(b=None, S='a'))               
Observable('obsCai', Ca(b=None, S='i'))
Observable('obsCaa', Ca(b=None, S='a'))
#Observable('obsE2Fu', E2F(b=None, S='u'))
#Observable('obsE2Fp', E2F(b=None, S='p'))
Observable('obspRBu', pRB(b=None, S='u'))
Observable('obspRBp', pRB(b=None, S='p'))
Observable('obspRBpp', pRB(b=None, S='pp'))
Observable('obsCdk2i', Cdk2(b=None, c=None, S='i'))
Observable('obsCdk2a', Cdk2(b=None, c=None, S='a'))
Observable('obsCa_Cdk2i', Ca(b=1, S='a') % Cdk2(b=1, S ='i'))
Observable('obsCa_Cdk2a', Ca(b=1, S='a') % Cdk2(b=1, S ='a'))
Observable('obsCa_Cdk2a_p27', Ca(b=1, S='i') % Cdk2(b=1, c=2, S ='i') % p27_p21(b=2, S='u'))
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
Observable('obsCbi', Cb(b=None, S='i'))
Observable('obsCba', Cb(b=None, S='a'))

#Initial Conditions

Parameter('Ca_0', 4.0e-05)           #Y2 in Iwamoto
Initial(Ca(b=None, S='a'), Ca_0)

Parameter('Cdk2_0', 1.35e01)           #Y4 in Iwamoto
Initial(Cdk2(b=None, c=None, S='i'),Cdk2_0)   

Parameter('Ca_Cdk2i_0', 4.0e-04)                       #Y8 in Iwamoto- unphosphorylated is inactive
Initial(Ca(b=1, S='a') % Cdk2(b=1, c=None, S='i'), Ca_Cdk2i_0)   

Parameter('Ca_Cdk2a_0', 1.0e-04)                       #Y9 in Iwamoto- phosphorylated is active
Initial(Ca(b=1, S='a') % Cdk2(b=1, c=None, S='a'), Ca_Cdk2a_0)   

Parameter('p27_21_0', 1.0e01)                   #Y10 in Iwamoto
Initial(p27_p21(b=None, S='u'),p27_21_0)

Parameter('Ca_Cdk2_p27_p21_0', 1.0e-04)                   #Y13 in Iwamoto
Initial(Ca(b=1, S='a') % Cdk2(b=1, c=2, S='a') % p27_p21(b=2, S='u'), Ca_Cdk2_p27_p21_0)

Parameter('CDC25_0', 1.0e01)				#THIS PARAMETER IS NOT CORRECT- ESTIMATED FOR NOW- NOT IN IWAMOTO PAPER	
Initial(CDC25(b=None, S='i'), CDC25_0)     

#Declare Rules

#E2F promotes the synthesis of Cyclin A 
#I can show synthesis, as shown below, but how do I show the promotion part????

#Cyclin A Synthesis and Degradation
synthesize(Ca(b=None, S = 'a'), kSCa)
degrade(Ca(b=None, S='a'),kDCa)

#Cyclin A binds to Cdk2 
Rule('Ca_Cdk2_bind',Ca(b=None, S='a') + Cdk2(b=None, c=None, S='i') <> Ca(b=1, S='a') % Cdk2(b=1, c=None, S ='i'), *[kcdk2f, kcdk2r])

#The Ca/Cdk2 complex equilibrates from inactive to active
equilibrate(Ca(b=1, S='a') % Cdk2(b=1, c=None, S='i'),Ca(b=1, S='a') % Cdk2(b=1, c=None, S='a'),[kfcacdk,krcacdk])

#The Ca/Cdk2 complex binds to p27_p21 - this complex is inactive
Rule ('Ca_Cdk2_p27_p21_bind',p27_p21(b=None, S='u') + Ca(b=1, S='a') % Cdk2(b=1, c=None, S='a') <> Ca(b=1, S='i') % Cdk2(b=1, c=2, S ='i') % p27_p21(b=2, S='u') , *[kfcp27_p21, krcp27_p21])

#Equilibrate CDC25
equilibrate(CDC25(b=None,S='i'),CDC25(b=None, S='a'),[kfCDC25, krCDC25])

#Equilibrate Cdh1
equilibrate(Cdh1(b=None,S='i'),Cdh1(b=None, S='a'),[kfCdh1, krCdh1])

#Previously, I thought synthesis should be shown as an enzyme reaction because it is promoted by a certain factor:
#catalyze(E2F(S='u'), 'b', Ca(S='i'), 'b', Ca(S='a'), (kfca, krca, kcca))  # calls the function catalyze to simulate activation of Cyclin E

#Parameters that bind to Ca
#bind_table([[ Ca(S= 'a')]
#	     [pRB(S = 'u'), pRB_rates, None]              #pRB inhibits Cyclin A
#	     [pRB(S = 'p'), pRBp_rates, None]             #pRBp inhibits Cyclin A 
#	     [cdk2(S = 'i'), None, cdk2_rates]])	  #cdk2 binds to Cylcin A to form a complex. 