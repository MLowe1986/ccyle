import TestG2M as m
import pysb.bng
from pysb.integrate import odesolve    #import odesolve to solve equations
import pysb.bng
import pylab as pl                     #import pylab to do plotting

t = pl.linspace(0, 200000)     #defines a time period (0 to 2000 seconds)
yout = odesolve(m.model, t)   #solves equations for the model over the specified time period


#Graph 1 
pl.ion()              #starts interactive mode
pl.figure()	      #creates a figure
#Say what to plot with a label for each
pl.plot(t,yout['obsCba'], label = "Active Cyclin B")
pl.plot(t,yout['obsCdk1i'], label = "Inactive Cdk1")
pl.plot(t,yout['obsCdk1a'], label = "Active Cdk1")
pl.plot(t,yout['obsCb_Cdk1i'], label = "Inactive Cyclin B/Cdk1 Complex")
pl.plot(t,yout['obsCb_Cdk1a'], label = "Active Cyclin B/Cdk1 Complex")
pl.plot(t,yout['obsCb_Cdk1a_p27'], label = "Inactive Cyclin B/Cdk1/p27_p21 Complex")



#Create a legend for the graph
pl.legend()
pl.xlabel("Time (s)")
pl.ylabel("Molecules/cell")

#Display graph
pl.show()

t = pl.linspace(0, 200000)     #defines a time period (0 to 2000 seconds)
yout = odesolve(m.model, t)   #solves equations for the model over the specified time period

#Graph 2
#pl.ion()
#pl.figure()

#What to plot:
#pl.plot(t,yout['obsCDC25i'], label = "Inactive CDC25")
#pl.plot(t,yout['obsCDC25a'], label= "Active CDC25")

pl.show()