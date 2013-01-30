
import Test as m

from pysb.integrate import odesolve    #import odesolve to solve equations
import pylab as pl                     #import pylab to do plotting

t = pl.linspace(0, 2000000)     #defines a time period (0 to 20,000 seconds)
yout = odesolve(m.model, t)   #solves equations for the model over the specified time period

pl.ion()              #starts interactive mode
pl.figure()	      #creates a figure
#Say what to plot with a label for each
pl.plot(t,yout['obsAP1i'], label = "Inactive AP1")
pl.plot(t,yout['obsAP1a'], label = "Active AP1")
pl.plot(t,yout['obsGF'], label = "Growth Factor")

#Create a legend for the graph
pl.legend()
pl.xlabel("Time (s)")
pl.ylabel("Molecules/cell")

#Display graph
pl.show()