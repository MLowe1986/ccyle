
from pysb import *
from pysb.macros import catalyze, degrade, equilibrate, bind, bind_table
    
# instantiate a model
Model()

#Declare Monomers
Monomer('Skp2', ['b','S'], {'S':['a', 'd']}) #Skp2 active form and degraded place holder
Monomer('Cdh1', ['b', 'S'], {'S':['i', 'a']})#Cdh1 active form and incative form
Monomer('pRB', ['b', 'S'],{'S':['u', 'p', 'pp']})  # Retinoblastoma protein unphosphorylated, phosphorylated once, and phosphorylated twice. 
Monomer('Ce',['b', 'd', 'S'], {'S':['i', 'a', 'd']})  #Cyclin E, inactive, active, and degraded (a place holder to indicate complete degradation of CyclinE)
Monomer('Cdk2',['b', 'c', 'S'], {'S':['i', 'a']})
Monomer('p27_p21',['b', 'S'], {'S':['u', 'p']})
Monomer('Wee1',['b', 'S'], {'S':['a', 'i']})
Monomer('CDC25', ['b', 'S'], {'S':['a', 'i']})
Monomer('E2F', ['b', 'S'], {'S':['a', 'i']})


#DECLARE PARAMETERS: 


# STEP 1: Parameters for Activation of Cyclin E by E2F 
Parameter('kSCe', 7.5e-02)     #This is k5 in Iwamoto paper


#STEP 2: Parameters for Inhibition of Cyclin E activation by pRB and pRBp. Calculated from Ki9 and Ki10 in Gerard. 
Parameter('kfpRB', 1.66e-06)
Parameter('krpRB', 1.00e-01)
Parameter('kfpRBp', 8.30e-08)
Parameter('krpRBp', 1.00e-01)


#STEP 3: Parameters for Degradation of Cyclin E
  #Iwamoto
Parameter('kDCe', 2.5e-03)    #This is k6 in Iwamoto paper
  #Gerard
Parameter('kfCe', 9.1e-08)    #(s^-1pL/mol) - Calculated from Km of 2microMol in Gerard and assuming krCe = 10^-1s^-1 and kcCe = 10^-2s^-1          
Parameter('krCe', 1e-1)       #Estimated
Parameter('kcCe', 1e-2)       #Estimated


#STEP 4: Parameters for Degradation of Skp2 
Parameter('kfcdh1', 4.57e-07)
Parameter('krcdh1', 1.00e-01)
Parameter('kccdh1', 1.00e-02) 



# STEP 5: Cyclin E/Cdk2 binding parameters: 
Parameter('kcdk2f',1.25e-03)
Parameter('kcdk2r', 2.5e-04)  #forward and reverse rates respectively for formation of CyclinE/Cdk2 complex (k7 and k8 from Iwamoto- CHECK THIS BY COMPARING TO EQUATIONS IN PAPER!!)

#STEP 6: Cdk2 irreversible flux parameter)
Parameter('kCdk2', 5.0e-04) # k16 in Iwamoto


#STEP 7: 

#CyclinE/Cdk2 Complex Equilibrate Parameters:
###I AM NOT SURE WHAT THESE ARE IF THEY ARE CORRECT
Parameter('kfcecdk', 2.5e-02)     #k22 from Iwamoto = NOT SURE
Parameter('krcecdk', 1.75e-03)    #k23 from Iwamoto = NOT SURE

#Parameters for self-catalysis of Ce/Cdk2(active) formation:
Parameter('kfcecom',    )   #WHERE ARE THESE PARAMETERS IN IWAMOTO PAPER?????
Parameter('krcecom',    )    
Parameter('kccecom',     )


Parameter('kfcp27_p21', 1.0e-02)     #k35 from Iwamoto
Parameter('krcp27_p21', 1.75e-04)    #k25 from Iwamoto

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
Observable('obsCe_Cdk2a_p27', Ce(b=1, S='i') % Cdk2(b=1, c=2, S ='i') % p27_p21(b=2, S='u'))
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

Parameter('Ce_Cdk2i_0', 1.0e-03)                       #Y6 in Iwamoto
Initial(Ce(b=1, S='a') % Cdk2(b=1, c=None, S='i'), Ce_Cdk2i_0)   

Parameter('Ce_Cdk2a_0', 1.0e-03)                       #Y7 in Iwamoto
Initial(Ce(b=1, S='a') % Cdk2(b=1, c=None, S='a'), Ce_Cdk2a_0)   

Parameter('p27_21_0', 1.0e01)                   #Y10 in Iwamoto
Initial(p27_p21(b=None, S='u'),p27_21_0)

Parameter('Ce_Cdk2_p27_p21_0', 1.0e00)                   #Y12 in Iwamoto
Initial(Ce(b=1, S='i') % Cdk2(b=1, c=2, S='i') % p27_p21(b=2, S='u'), Ce_Cdk2_p27_p21_0)

Parameter('CDC25_0', 1.0e01)				#THIS PARAMETER IS NOT CORRECT- ESTIMATED FOR NOW- NOT IN IWAMOTO PAPER	
Initial(CDC25(b=None, S='i'), CDC25_0)     

#DECLARE RULES: 



# STEP 1: E2F promotes the activation of Cyclin E  
    # FIXME: (Select a model)
    #Iwamoto Model: 
Rule('E2F_Ce_trans', E2F() >> Ce, kSCe) #Iwamoto models activation of Cyclin E as a transformation from E2F to Cyclin E
    #Gerard Model: 
#catalyze(E2F(S='a'), 'b', Ce(S='i'), 'b', Ce(S='a'), (kfCe, krCe, kcCe))  #Cyclin E Activation is promoted by E2F, which can be modeled as a catalysis reaction. Parameters from Gerard also seem to show this as a transformation...(???)



#STEP 2: pRB and pRBp inhibit Cyclin E activation by binding to inactive Cyclin E
   #Iwamoto Model: Does not include this inhibition step
   #Gerard Model:
#Parameters that bind to Ce
bind_table([[Ce(S= 'a')]
	     [pRB(S = 'u'), (kfpRB, krpRB)]                         #pRB inhibits Cyclin E   
	     [pRB(S = 'p'), (kfprBp, krpRBp)]], 'b', 'b')            #pRBp inhibits Cyclin E  


# STEP 3: Skp2 promotes the degradation of Cyclin E: 
     #FIXME: (Select a model)
     #Iwamoto Model:
#degrade(Ce(b=None, S='a'),kDCe)
     #Gerard Model:
catalyze(Skp2(S = 'a'), 'b', Ce(S='a'), 'b', Ce(S='d'), (kfce, krce, kcce))  # calls the function catalyze to simulate activation of Cyclin E


#STEP 4: Cdh1 promotes the degradation of Skp2
    #Iwamoto Model:Does not include this degradation step
    #Gerard Model: 
catalyze(Cdh1(S= 'a'), 'b', Skp2(S='a'), 'b', Skp2(S='d'), (kfcdh1, krcdh1, kccdh1)) 


#STEP 5:Cyclin E binds to Cdk2 
   #Iwamoto Model: Same as Gerard model
   #Gerard Model:
Rule('Ce_Cdk2_bind',Ce(b=None, S='a') + Cdk2(b=None, c=None, S='i') <> Ce(b=1, S='a') % Cdk2(b=1, c=None, S ='i'), *[kcdk2f, kcdk2r])



#STEP 6: Cdk2 leave the Ce/Cdk2 complex via irreversible flux
  #FIXME: IS THIS THE CORRECT EQUATION FOR THIS STEP????
  #Iwamoto Model: 
  Rule('Ce_Cdk2_trans', Ce(b=1, S='a') % Cdk2(b=1, c=None, S ='i') >> Cdk2(b=None, c=None, S='i'), kCdk2) 
  #Gerard Model: This step is not shown
  
  
#STEP 7: Ce/Cdk2 becomes active (phosphorylated)
  #FIXME: Select a model
   #Iwamoto: #The active (phosphorylated) form of the Ce/Cdk2 complex self activates
catalyze(Ce(b=1, d=None, S='a') % Cdk2(b=1, c=None, S='a'), 'd', Ce(b=1, d=None S='a') % Cdk2(b=1, c=None, S='i'), 'd', Ce(b=1, d=None, S='a') % Cdk2(b=1, c=None, S='a'), (kfcecom, krcecom, kccecom))
   #Gerard: The Ce/Cdk2 complex equilibrates from inactive to active 
equilibrate(Ce(b=1, S='a') % Cdk2(b=1, c=None, S='i'),Ce(b=1, S='a') % Cdk2(b=1, c=None, S='a'),[kfcecdk,krcecdk]) #currently these parameters are Iwamoto parameters...


#STEP 8: CDC25 equilibrates from active to inactive. This reaction is catalyzed by the active Ce/Cdk2 complex. 
  #Iwamoto: This step is not shown
  # Gerard: 
catalyze(Ce(b=1, d=None, S='a') % Cdk2(b=1, c=None, S='a'), 'd', CDC25(b=None, S='a'), 'b', CDC25(b=None, S='a'), (kfcdc, krcdc, kccdc))
equilibrate(CDC25(b=None,S='i'),CDC25(b=None, S='a'),[kfCDC25, krCDC25])



#The Ce/Cdk2 complex binds to p27_p21 - this complex is inactive
Rule ('Ce_Cdk2_p27_p21_bind',p27_p21(b=None, S='u') + Ce(b=1, S='a') % Cdk2(b=1, c=None, S='a') <> Ce(b=1, S='i') % Cdk2(b=1, c=2, S ='i') % p27_p21(b=2, S='u') , *[kfcp27_p21, krcp27_p21])

#Phosphorylate p27_21
equilibrate(p27_p21(b=None,S='u'),p27_p21(b=None, S='p'),[kfp27_p21, krp27_p21])



