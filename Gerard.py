
# import the pysb module and all its methods and functions
from pysb import *
from pysb.macros import catalyze, bind

def catalyze(enz, sub, site, state1, state2, kf, kr, kc):   # function call
    """2-step catalytic process"""                          # reaction name
    r1_name = '%s_assoc_%s' % (enz.name, sub.name)           # name of association reaction for rule
    r2_name = '%s_diss_%s' % (enz.name, sub.name)           # name of dissociation reaction for rule
    E = enz(b=None)                                         # define enzyme state in function
    S = sub({'b': None, site: state1})                      # define substrate state in function
    ES = enz(b=1) % sub({'b': 1, site: state1})             # define state of enzyme:substrate complex
    P = sub({'b': None, site: state2})                      # define state of product
    Rule(r1_name, E + S <> ES, kf, kr)                      # rule for enzyme + substrate association (bidirectional)
    Rule(r2_name, ES >> E + P, kc)                          # rule for enzyme:substrate dissociation  (unidirectional)
    
    
# instantiate a model
Model()

# declare monomers

# ???do we enter in factors that are going to promote synthesis etc, but not necessarily bind as monomers? 

Monomer('GF',['b'])        #Growth Factor that activates AP1
Monomer('AP1', ['b'])      #Transcription factor that promotes the synthesis of cyclin D.
Monomer('Ca', ['b'])       #Cyclin A
Monomer('Cb', ['b'])       #Cyclin B
Monomer('Cd', ['b'])       #Cyclin D
Monomer('Ce', ['b'])       #Cylin E
Monomer('Cdc20', ['b', 'S'], {'S':['a', 'i']})   #active and inactive forms of Cdc20
Monomer('Cdh1', ['b', 'S'], {'S':['a', 'i']})    #active and inactive forms of Cdh1
Monomer('E2F', ['b', 'S'], {'S': ['t', 'p']})   #transcription factor E2F and the phosphorylated form of E2F

# ???For the capital letters- these are not Monomers, so don't declare as such. But how do we declare an active and inactive state for a complex? 
 
 #A.  
  
Monomer('Ma', ['b', 'S'], {'S':['a', 'i']})  # ??For these ones... Ma is for cyclin A/Cdk2 complex...        
                                              #don't we need to define each individual and use the rule to form the complex? Then the monomers that we define will be different than what is in their paper.
   #OR 
     #a. 
Monomer('Cdk2', ['b'])  #with the rule below 
 
# B. 
Monomer('Mb', ['b', 'S'], {'S':['a', 'i']})   #For these ones... Mb is for cyclin B/Cdk1
  
#OR 
     #b. 
Monomer('Cdk1', ['b'])  #with the rule below 


# C. 
Monomer('Md', ['b', 'S'], {'S':['a', 'i']})   # For these ones... Md is for cyclin D/Cdk4-6

#OR 

    #c. 
Monomer('Cdk4_6', ['b'])  #with the rule below 
    
#D.
Monomer('Me', ['b', 'S'], {'S':['a', 'i']})  # For these ones... Me is for cyclin E/Cdk2

#OR 

    #d.
Monomer('Cdk2', ['b'])  #with the rule below   * This is the same as a above, so we don't need to declare the monomer twice. 
  
#E  Map27, Mdp27, Mep27, Mbp27  - These are defined in paper as a variable, but they are not monomers- do we define p27 instead and then create rules to join to the complexes? 
  
#e.  
Monomer('p27', ['b', 'S'], {'S':['u', 'p']})   #with rules below  Cyclin-dependent kinase inhibitor p27 (or p21), unphosphorylated or phosphorylated
  
  
# ?For phosphatase Cdc25, the paper creates multiple variables of the same monomer, depending on which monomer it acts above (and also if it is active or inactive. Eg: Pa acts on cyclinA/Cdk2, Pb acts on cyclinB/Cdk1. Can we just declare 1 monomer, but many rules? Pe acts on cyclinE/Cdk2

Monomer('P', ['b'], {'S':['a', 'i']}   # an active, phosphorylated form and an inactive, dephosphorylated form of phosphase Cdc25, which acts on cyclinA/Cdk2, cyclinB/Cdk1, and cyclinE/Cdk2. 

Monomer('pRB', ['b', 'S'],{'S':['u', 'p', 'pp']})    #can it have more than two states? Retinoblastoma protein unphosphorylated, phosphorylated once, and phosphorylated twice. 

#pRBcl and pRBc2 are listed as variables in the paper, but these are just complexes of pRB/E2F and pRBp/E2F defined by the rules below. 

Monomer('Skp2', ['b'])#F-box protein Skp2

Monomer('Wee1', ['b'], {'S':['u', 'p']})    #Active, dephosphorylated form of kinase Wee1 and inactive, phosphorylated form of kinase Wee1

Monomer('ATR', ['b'])    # Kinase ATR
Monomer('Cdc45', ['b'])   #DNA anchor factor
Monomer('Chk1', ['b'])    #Kinase Chk1
Monomer('Pol', ['b'])    #DNA polymerase alpha
Monomer('Primer', ['b'])  #RNA primer
Monomer('Bn', ['b'])    #Nuclear form of CLOCK-BMAL1 complex
Monomer('Mw', ['b'])    #Wee1 mRNA induced by the circadian clock 

# input the parameter values
#?? Are we going to use the same units? 
Parameter('kf', 1.0e-07)
Parameter('kr', 1.0e-03)
Parameter('kc', 1.0)

# now input the rules

#COULD WE SIMPLIFY THE RULES BY ADDING IN A BIND TABLE? 

Rule('cplx_bind', 
     Ca(b=1) % Cdk2(b=1, bp=None) + P27(b=None) <>
     Ca(b=1) % Cdk2(b=1, bp=2) % P27(b=2),
     KF, KR)

bind(Ca(b=1) % Cdk2(b=1), 'bp' P27(), 'b', [KF, KR])

A. Rule('Ca_Cdk2_bind', Ca(b=None) + Cdk2(b=None) <> Ca(b=1) % Cdk2(b=1), *[kcom3, kdecom3]) *** HOW DO WE DEFINE A STATE FOR THE COMPLEX- ACTIVE AND INACTIVE? OR DO WE MAKE ANOTHER RULE WHERE THE COMPEXT IS ACTIVATED? 

B. Rule('Cb_Cdk1_bind', Cb(b=None) + Cdk1(b=None) <> Cb(b=1) % Cdk1(b=1), *[kcom4, kdecom4])

C. Rule('Cd_Cdk4_6_bind', Cd(b=None) + Cdk4_6(b=None) <> Cd(b=1) % Cdk4_6(b=1), *[kcom1, kdecom1])

D. Rule('Ce_Cdk2bind', Ce(b=None) + Cdk2(b=None) <> Ce(b=1) % Cdk2(b=1), *[kcom2, kdecom2])

E. Rule('Ca_Cdk2_p27_bind', Ca(b=1) % Cdk2(b=1) + p27(b=None, S='u') <> Cd(b=2) % Cdk2(b=2) % p27(b=2, S='u'), *[kc5, kc6])

Rule('Cd_Cdk4_6_P27_bind', Cd(b=1) % Cdk4_6(b=1) + p27(b=None, S='u') <> Cd(b=2) % Cdk4_6(b=2) % p27(b=2, S='u'), *[kc1, kc2])

Rule('Ce_Cdk2_P27_bind', Ce(b=1) % Cdk2(b=1) + p27(b=None, S='u') <> Ce(b=2) % Cdk2(b=2) % p27(b=2, S='u'), *[kc3, kc4]

Rule('Cb_Cdk1_P27_bind', Cb(b=1) % Cdk1(b=1) + p27(b=None, S='u') <> Cb(b=2) % Cdk1(b=2) % p27(b=2, S='u'), *[kc7, kc8])

Rule('pRB_E2F__bind', pRB(b=None, S= 'u') + E2F(b=None, S = 't') <> pRB(b=1, S='u') % E2F(b=1, S= 't'), *[kpc1, kpc2])  #??Does the phosphorylated or unphosphorylated form react with pRB? 

Rule('pRBp_E2F__bind', pRB(b=None, S= 'p') + E2F(b=None, S = 't') <> pRB(b=1, S='p') % E2F(b=1, S= 't'), *[kpc3, kpc4])  #??Does the phosphorylated or unphosphorylated form react with pRB? 


#** HOW WOULD WE SHOW GROWTH FACTOR PROMOTING SYNTHESIS??
catalyze(GF, AP1, 'S', 'u', 't', kf, kr, kc)



# initial conditions
Parameter('C8_0', 1000)
Parameter('Bid_0', 10000)
Initial(C8(b=None), C8_0)
Initial(Bid(b=None, S='u'), Bid_0)

# Observables
Observable('obsC8', C8(b=None))
Observable('obsBid', Bid(b=None, S='u'))
Observable('obstBid',Bid(b=None, S='t'))v
