"""
Python script to test the qutip python library.

I define the steady state solutions of the Bloch equations for
a magnetization M = [Mx, My, Mz]. 

I should be able to reproduce these curves using the steady state
solver of the qutip library.

"""

from __future__ import division
import qutip as qt
import numpy as np
import matplotlib.pyplot as plt

# Setting parameters for the simulation

w0 = 2 * np.pi
w1 = 0.5 * 2 * np.pi
T1 = 0.5
T2 = 0.5
M0 = 1

# Definition of the steady-state solutions of the Bloch equations for magnetization M = [Mx, My, Mz]
# w is the frequency of the applied RF-field

def Mx(w):
    nominator = w1 * T2**2 * (w0 - w)
    denominator = 1 + w1**2 * T1 * T2 + (T2 * (w0 - w))**2

    return nominator/denominator * M0

def My(w):
    nominator = w1 * T2
    denominator = 1 + w1**2 * T1 * T2 + (T2*(w0 - w))**2
    
    return nominator/denominator * M0

def Mz(w):
    nominator = 1 + T2**2 * (w0 - w)**2
    denominator = 1 + w1**2 * T1 * T2 + (T2 * (w0 - w))**2

    return nominator/denominator * M0



# Calculating steady-state solutions using the qutip library

def H(w):
    """ Hamiltonian of a spin 1/2 in B0 + Brf using the RWA. """
    return (w0 - w)/2 * qt.sigmaz() + w1/2 * qt.sigmax()

# Collsapse operators which should describe the relaxation of the system
c1 = np.sqrt(1.0/T1) * qt.destroy(2)
c2 = np.sqrt(1.0/T2)/2 * qt.destroy(2)
c3 = np.sqrt(1.0/(T1*T2)) * qt.destroy(2)
c4 = np.sqrt(1.0/T2) /2 * qt.create(2)
c5 = np.sqrt(1.0/(T1*T2)) * qt.create(2)


wHVals = np.linspace(-2 * np.pi, 6 * np.pi, 31)
MxH = [qt.expect(qt.sigmax(), qt.steadystate(H(w), [c1,c2,c3], maxiter=10)) for w in wHVals]
MyH = [qt.expect(qt.sigmay(), qt.steadystate(H(w), [c2,c5], maxiter=10)) for w in wHVals]
MzH = [qt.expect(qt.sigmaz(), qt.steadystate(H(w), [c1], maxiter=10)) for w in wHVals]

wVals = np.linspace(-2 * np.pi, 6 * np.pi, 100)

# Plot results of the algebraic formula and the solution given by qutip

plt.ylim([-1.0, 1.0])

plt.plot(wVals, Mx(wVals), label='Mx')
plt.plot(wVals, My(wVals), label='My')
plt.plot(wVals, Mz(wVals), label='Mz')

plt.plot(wHVals, MxH, 'bo', label='MxH')
plt.plot(wHVals, MyH, 'go', label='MyH')
plt.plot(wHVals, MzH, 'ro', label='MzH')

plt.legend(loc=3, numpoints=1)
plt.show()


"""
psi0 = qt.basis(2, 0)
tlist = np.linspace(0, 15, 100)

result = qt.mesolve(H, psi0, tlist, c_op_list, [])

expectZ = qt.expect(qt.sigmaz(), result.states)
expectX = qt.expect(qt.sigmax(), result.states)
expectY = qt.expect(qt.sigmay(), result.states)

plt.plot(tlist, expectZ, label=r'$<\sigma_Z >$')
plt.plot(tlist, expectX, label=r'$<\sigma_X >$')
plt.plot(tlist, expectY, label=r'$<\sigma_Y >$')
plt.legend()
plt.show()
"""

