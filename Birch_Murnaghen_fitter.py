#!/usr/bin/env python

import csv

from pylab import * #this includes numpy as np!
from scipy.optimize import leastsq

csv.register_dialect('ssv', delimiter=' ', skipinitialspace=True)

# pulls the data from input.txt
data = []
with open('input.txt', 'r') as f:
    reader = csv.reader(f, 'ssv')
    for row in reader:
        floats = [float(column) for column in row]
        data.append(floats)
fullData = np.array(data)

# Assigns first column to v, and second column to e, respectively
v = fullData[:,0]
e = fullData[:,1]

# Make a vector to evaluate fits on with a lot of points so it looks smooth
vfit = np.linspace(min(v),max(v),100)

# Fit a parabola to the data
# y = ax^2 + bx + c
a,b,c = polyfit(v,e,2) #this is from pylab

# Initial guesses.
v0 = -b/(2*a)
e0 = a*v0**2 + b*v0 + c
b0 = 2*a*v0
bP = 2

# Create the equation of state function
def BirchMurnaghan(parameters,vol):

    E0 = parameters[0]
    B0 = parameters[1]
    BP = parameters[2]
    V0 = parameters[3]

    E = E0 + (((9*V0*B0)/16)*(((((V0/vol)**(2/3))-1)**3)*BP + ((((V0/vol)**(0.666))-1)**2)*(6-(4*((V0/vol)**(2/3))))))
    return E

# Define an objective function that will be minimized
def objective(pars,y,x):
    #we will minimize this function
    err =  y - BirchMurnaghan(pars,x)
    return err

x0 = [e0, b0, bP, v0] #initial guesses in the same order used in the Murnaghan function

murnpars, ier = leastsq(objective, x0, args=(e,v)) #this is from scipy

# Create a figure summarizing the results
plot(v,e,'ro')
plot(vfit, a*vfit**2 + b*vfit + c,'--',label='parabolic fit')
plot(vfit, BirchMurnaghan(murnpars,vfit), label='Birch-Murnaghan fit')
xlabel('Volume ($\AA^3$)')
ylabel('Ground state energy (eV)')
legend(loc='best')

# Add text to the figure in figure coordinates
ax = gca()
# text(0.4,0.7,'V0 = %1.9f $\AA^3$' % murnpars[3],
#      transform = ax.transAxes)
# text(0.4,0.5,'E0 = %1.9f $\AA^3$' % murnpars[0],
#      transform = ax.transAxes)
# text(0.4,0.4,'BP = %1.9f $\AA^3$' % murnpars[2],
#      transform = ax.transAxes)
# text(0.4,0.6,'B0 = %1.9f eV/$\AA^3$' % murnpars[1], transform = ax.transAxes)
text(0.35,0.6,'Bulk modulus = %1.2f GPa' % (murnpars[1]*160.21773), transform = ax.transAxes)
savefig('output.png')
show()


print 'initial guesses  : ',x0
print 'fitted parameters: ', murnpars
