# import the pysb module and all its methods and functions
from pysb import *
from pysb.macros import catalyze, degrade, equilibrate, bind
    
# instantiate a model
Model()


#Declare Monomers
Monomer('GF',['b'])                               #Growth Factor that activates synthesis of AP1
Monomer('AP1', ['b', 'S'],{'S':['i','a']})       #Transcription factor that promotes the synthesis of cyclin D.
Monomer('Cd', ['b', 'S'], {'S':['i', 'a']})       #Cyclin D
Monomer('E2F', ['b', 'S'], {'S': ['u', 'p']})     #transcription factor E2F and the phosphorylated form of E2F
Monomer('pRB', ['b', 'S'],{'S':['u', 'p', 'pp']})  # Retinoblastoma protein unphosphorylated, phosphorylated once, and phosphorylated twice. 
Monomer('Cdk4_6', ['b'], {'S':['a', 'i']})          #Cdk4_6 are treated as one monomer
Monomer('p27_p21',['b'])
Monomer('Ce',['b'])
Monomer('Cdk2',['b'], {'S':['a', 'i']})
Monomer('Ca',['b'])

#Paramaters taken from Iwamoto paper:

Parameter(







#AP1 Activation and Degradation Parameters 
Parameter('kfap1', )	   #AP1 activation forward rate 
Parameter('krap1', )       #AP1 activation reverse rate 
Parameter('kcap1', )       #AP1 activation catalytic rate
Parameter('kdap1', )       #AP1 Degradation Rate

#CyclinD Activation by AP1 Parameters
Parameter('kfcdA', )	#Cyclin D activation by AP1 forward rate 
Parameter('krcdA', )    #CyclinD activation by AP1 reverse rate 
Parameter('kccdA', )    #CyclinD activation by AP1 catalytic rate 

#Cyclin D Activation Parameters by E2F
Parameter('kfcdE', )	   #Cyclin D activation by E2F forward rate 
Parameter('krcdE', )       #CyclinD activation by E2F reverse rate 
Parameter('kccdE', )       #CyclinD activation by E2F catalytic rate 

Parameter('kdddA', )    #CyclinD Degradation Rate

#forward and reverse rates for CyclinD inhibition by pRB 
pRBCd_rates = [   ,   ]
#forward and reverse rates for CyclinD inhibition by pRBp 
pRBpCd_rates = [   ,   ]

#forward and reverse rates for CyclinD/Cdk4_6 complex formation  
Cdk4_6_rates = [   ,   ]
 
#forward and reverse rates for CyclinD/Cdk4_6 complex equlibration
Parameter('kfcdcdk', )
Parameter('krcdcdk', ) 

#forward and reverse rates for CyclinD/Cdk4_6/p27_p21 complex formation
Parameter('kfp27',   )
Parameter('krp27',   ) 

#forward and reverse rates for phosphorylation of pRB
Parameter('kfprb',   )
Parameter('krprb',   )
 
#forward and reverse rates for phosphorylation of pRBp
Parameter('kfprbp',  )
Parameter('krprbp',  ) 

#forward and reverse rates for phosphorylation of E2F
Parameter('kfEp',   )
Parameter('krEp',   )


 
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
catalyze(AP1(S='a'), 'b', Cd(S='i'), 'b', Cd(S='a'), (kfcdA, krcdA, kccdA))  # calls the function catalyze to simulate synthesis of AP1 
  
#E2F activates Cyclin D  
catalyze(E2F(S='u'), 'b', Cd(S='i'), 'b', Cd(S='a'), (kfcdE, krcdE, kccdE))  # calls the function catalyze to simulate synthesis of AP1   

#Bind table: cdk4_6 binds to Cd to form a complex, which will then become activated in an equilibration step. pRB and pRBp bind to Cd to inhibit it. 
bind_table_Cd([[ Cd(S = 'i'), Cd(S = 'a')]
	     [pRB(S = 'u'), pRBCd_rates, None]              #pRB inhibits Cyclin D 
	     [pRB(S = 'p'), pRBpCd_rates, None]             #pRBp inhibits Cyclin D 
	     [Cdk4_6(S = 'i'), None, Cdk4_6_rates]])      #cdk4_6 and Cyclin D form a complex
	     
# The Cd/Cdk4_6	complex equilibrates from inactive to active 
equilibrate(Cd(b=1)%Cdk4_6(b=None,S='i'),Cd(b=1)%Cdk4_6(b=None,S='a'),[kfcdcdk,krcdcdk])

#The Cd/Cdk4_6 complex binds to p27_p21
Rule ('Cd_Cdk4_6_p27_p21_bind',p27_21(b=None) + Cd(b=1) % Cdk4_6(b=1, S='a') <> Cd(b=1) % Cdk4_6(b=1,b=2, S ='a') % p27_p21(b=2) , *[kfp27, krp27])

#pRB is phosphorylated to become pRBp
equilibrate(pRB(b=None,S='u'),pRB(b=None,S='p'),[kfprb,krprb]

#pRBp is phosphorylated to become pRBpp
equilibrate(pRB(b=None,S='p'),pRB(b=None,S='pp'),[kfprbp,krprbp]

#E2F is phosphorylated to become E2Fp
equilibrate(E2F(b=None,S='u'),E2F(b=None, S='p'),[kfEp, krEp]

catalyze(Ce(b=1)%Cdk2(b=1,S='a'), 'b', pRB(S='p'), 'b', pRB(S='pp'), (kfCe, krCe, kcaCe))

catalyze(Ca(b=1)%Cdk2(b=1,S='a'), 'b', E2F(S='u'), 'b', E2F(S='p'), (kfCa, krCa, kcaCa))

catalyze(Cd(b=1)%Cdk4_6(b=1,S='a'), 'b', pRB(S='u'), 'b', pRB(S='p'), (kfCd, krCd, kcaCd))

catalyze(Cd(b=1) % Cdk4_6(b=1,b=2 S ='a') % p27_p21(b=2)), 'b', pRB(S='u'), 'b', pRB(S='p'), (kfCdp, krCdp, kcaCdp))
  
#Rule for degradation of AP1 
degrade(AP1(), kdap1)

#Degradation of CyclinD
degrade(Cd(),kdCd)


#Rules replaced above by Bind: 
#pRB inhibits Cyclin D 
#Rule('pRB_Cd__bind', pRB(b=None, S= 'u') + Cd(b=None, S = 'i') <> pRB(b=1, S='u') % Cd(b=1, S= 'i'), *[kfi7, kri7])  #this can use bind
#pRBp inhibits Cyclin D 
#Rule('pRBp_Cd__bind', pRB(b=None, S= 'p') + Cd(b=None, S = 'i') <> pRB(b=1, S='p') % Cd(b=1, S= 'i'), *[kfi8, kri8])  #this can use bind

#Rule('Cd_Cdk4_6_bind', Cd(b=None, S='a') + Cdk4_6(b=None, S='i') <> Cd(b=1) % Cdk4_6(b=1, S= 'i'), *[kcom1, kdecom1])   #this can use bind