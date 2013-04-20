
import TestG1S as m
import pysb.bng
from pysb.integrate import odesolve    #import odesolve to solve equations
import pysb.bng
import pylab as pl                     #import pylab to do plotting

t = pl.linspace(0, 2000)     #defines a time period (0 to 2000 seconds)
yout = odesolve(m.model, t)   #solves equations for the model over the specified time period


#Graph 1 
pl.ion()              #starts interactive mode
pl.figure()	      #creates a figure
#Say what to plot with a label for each
pl.plot(t,yout['obsCea'], label = "Active Cyclin E")
pl.plot(t,yout['obsCdk2i'], label = "Inactive Cdk2")
pl.plot(t,yout['obsCdk2a'], label = "Active Cdk2")
pl.plot(t,yout['obsCe_Cdk2i'], label = "Inactive Cyclin E/Cdk2 Complex")
pl.plot(t,yout['obsCe_Cdk2a'], label = "Active Cyclin E/Cdk2 Complex")
pl.plot(t,yout['obsCe_Cdk2a_p27'], label = "Active Cyclin E/Cdk2/p27_p21 Complex")
pl.plot(t,yout['obsp27_p21u'], label = "Unphosphorylated p27_p21")
pl.plot(t,yout['obsp27_p21p'], label = "Phosphorylated p27_p21")


#Create a legend for the graph
pl.legend()
pl.xlabel("Time (s)")
pl.ylabel("Molecules/cell")

#Display graph
pl.show()

t = pl.linspace(0, 200000)     #defines a time period (0 to 2000 seconds)
yout = odesolve(m.model, t)   #solves equations for the model over the specified time period

#Graph 2
pl.ion()
pl.figure()

#What to plot:
pl.plot(t,yout['obsCDC25i'], label = "Inactive CDC25")
pl.plot(t,yout['obsCDC25a'], label= "Active CDC25")

pl.show()