# import the pysb module and all its methods and functions
from pysb import *
from pysb.macros import catalyze, degrade, equilibrate, bind
    
# instantiate a model
Model()


#Declare Monomers
Monomer('Cdh1',['b', 'S'],{'S':['i', 'a']})                            
Monomer('Skp2', ['b'])
Monomer('E2F', ['b', 'S'], {'S': ['u', 'p']})     #transcription factor E2F and the phosphorylated form of E2F
Monomer('pRB', ['b', 'S'],{'S':['u', 'p', 'pp']})  # Retinoblastoma protein unphosphorylated, phosphorylated once, and phosphorylated twice. 
Monomer('Ce',['b', 'S'], {'S':['i', 'a']})
Monomer('Cdk2',['b', 'S'], {'S':['i', 'a']})
Monomer('p27_p21',['b', 'S'], {'S':['u', 'p']})
Monomer('Wee1',['b'], {'S':['a', 'i']})


#Declare Parameters

#Parameter for Synthesis of Cyclin E 
Parameter('kSCe', 7.5e-02)	    #This is k5 in Iwamoto paper

#Parameter for Degradation of Cyclin E
Parameter('kDCe', 2.5e-03)

#Cyclin E/Cdk2 binding parameters: 
cdk2_rates = [1.25e-03, 2.4e-04]  #forward and reverse rates respectively for formation of CyclinE/Cdk2 complex (k7 and k8 from Iwamoto- CHECK THIS BY COMPARING TO EQUATIONS IN PAPER!!)

#Cyclin E/Cdk2/p27_p21 binding parameters:
p27_rates = [2.5e-02, 1.75e-03]  #forward and reverse rates respectively for formation of CyclinE/Cdk2/p27_p21 complex (k22 and k23 from Iwamoto- CHECK THIS BY COMPARING TO EQUATIONS IN PAPER!!)

#Observables
Observable('obsCdh1i', Cdh1(b=None, S='i'))    #observe the inactive form of AP1
Observable('obsCdh1a', Cdh1(b=None, S='a'))    #observe the active form of AP1
Observable('obsSkp2', Skp2(b=None))             #observe Growth factor
Observable('obsCei', Ce(b=None, S='i'))
Observable('obsCea', Ce(b=None, S='a'))
Observable('obsE2Fu', E2F(b=None, S='u'))
Observable('obsE2Fp', E2F(b=None, S='p'))
Observable('obspRBu', pRB(b=None, S='u'))
Observable('obspRBp', pRB(b=None, S='p'))
Observable('obspRBpp', pRB(b=None, S='pp'))
Observable('obsCdk2i', Cdk2(b=None, S='i'))
Observable('obsCdk2a', Cdk2(b=None, S='a'))
Observable('obsp27_p21u', p27_p21(b=None, S='u'))
Observable('obsp27_p21p', p27_p21(b=None, S='p'))
Observable('obsWee1i', Wee1(b=None, S='i'))
Observable('obsWee1a', Wee1(b=None, S='a'))

#Declare Rules

#E2F promotes the synthesis of Cyclin E  
#I can show synthesis, as shown below, but how do I show the promotion part????
Rule('Synthesize_Ce', None >> Ce(), kSCe)


#Previously, I thought synthesis should be shown as an enzyme reaction because it is promoted by a certain factor:
#catalyze(E2F(S='u'), 'b', Ce(S='i'), 'b', Ce(S='a'), (kfce, krce, kcce))  # calls the function catalyze to simulate activation of Cyclin E

#Parameters that bind to Ce
bind_table_Ce([[ Ce(S= 'a')]
	     [pRB(S = 'u'), pRB_rates, None]              #pRB inhibits Cyclin E 
	     [pRB(S = 'p'), pRBp_rates, None]             #pRBp inhibits Cyclin E 
	     [cdk2(S = 'i'), None, cdk2_rates]])	  #cdk2 binds to Cylcin E to form a complex. 
	    
#The Ce/Cdk2 complex binds to p27_p21
Rule ('Ce_Cdk2_p27_p21_bind',p27_21(b=None) + Ce(b=1) % Cdke(b=1, S='a') <> Ce(b=1) % Cdk2(b=1,b=2, S ='a') % p27_p21(b=2) , *[p27_rates])


#The Ce/Cdk2 complex equilibrates from inactive to active
equilibrate(Ce(b=1)%Cdk2(b=None,S='i'),Ce(b=1)%Cdk2(b=None,S='a'),[kfcecdk,krcecdk])

equilibrate(p27_p21(b=None,S='u'),p27_p21(b=None, S='p'),[kfp27, krp27]

#Degradation of Cyclin E
#degrade(Ce(), kDCe)