import numpy as np
import scipy.constants

# This code is attributable to Jamal I. Mustafa
# in the unreleased software cmsPy version 0.0.
# Used with permission(?)

# energy conversions
Ha2eV = scipy.constants.value('Hartree energy in eV')
Ha2meV = Ha2eV * 1000
Ry2eV = scipy.constants.value('Rydberg constant times hc in eV')
Ha2Ry = 2.0
Ry2Ha = 0.5
Ha2K = scipy.constants.value('hartree-kelvin relationship')

# length conversions
Bohr2Angstrom = scipy.constants.value('Bohr radius') * 1e10

# epsilons
eps = np.finfo(np.float).eps
eps3 = 1e-3
eps6 = 1e-6
eps9 = 1e-9
eps12 = 1e-12
