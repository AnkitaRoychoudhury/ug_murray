# This file reads an SBML file using libSBML, 
#
# - expands all function definitions
# - expands all initial assignments
# - converts local parameter to global ones
# - then it goes ahead and write the ODE system
#   for use with scipy.integrate
# - it emits a function called simulateModel
#   that takes three parameters: t0, tend and numpoints
#   with that the model can be simulated as needed
# - finally the emitted function is called and the result plotted 
# - it is also written out into a file called generated.py 
#

from libsbml import *
import sys


def generateCodeForFile(filename, t0=0, tEnd=10, numPoints=1000):
  # 
  # read the SBML from file 
  # 
  doc = readSBMLFromFile(filename)
  if doc.getNumErrors(LIBSBML_SEV_FATAL):
    print('Encountered serious errors while reading file')
    print(doc.getErrorLog().toString())
    sys.exit(1)
    
  # clear errors
  doc.getErrorLog().clearLog()
  
  #
  # perform conversions
  #
  props = ConversionProperties()
  props.addOption("promoteLocalParameters", True)
  
  if doc.convert(props) != LIBSBML_OPERATION_SUCCESS: 
    print('The document could not be converted')
    print(doc.getErrorLog().toString())
    
  props = ConversionProperties()
  props.addOption("expandInitialAssignments", True)
  
  if doc.convert(props) != LIBSBML_OPERATION_SUCCESS: 
    print('The document could not be converted')
    print(doc.getErrorLog().toString())
    
  props = ConversionProperties()
  props.addOption("expandFunctionDefinitions", True)
  
  if doc.convert(props) != LIBSBML_OPERATION_SUCCESS: 
    print('The document could not be converted')
    print(doc.getErrorLog().toString())
    
  #
  # figure out which species are variable 
  #
  mod = doc.getModel()
  variables = {}
  
  for i in range(mod.getNumSpecies()): 
     species = mod.getSpecies(i)
     if species.getBoundaryCondition() == True or variables.has_key(species.getId()):
       continue
     variables[species.getId()] = []
  
  #
  # start generating the code by appending to bytearray
  #
  generated_code = bytearray('')
  generated_code.extend('from numpy import *\n')
  generated_code.extend('from matplotlib.pylab import *\n')
  generated_code.extend('from matplotlib.pyplot import *\n')
  generated_code.extend('from scipy.integrate import odeint \n')
  
  generated_code.extend('\n')
  generated_code.extend('def simulateModel(t0, tend, numPoints):\n')
  
  # write out compartment values 
  generated_code.extend('  \n')
  generated_code.extend('  #compartments\n')
  for i in range(mod.getNumCompartments()):
    element = mod.getCompartment(i)
    generated_code.extend('  {0} = {1}\n'.format(element.getId(), element.getSize()))
  
  # write out parameter values 
  generated_code.extend('  \n')
  generated_code.extend('  #global parameters\n')
  for i in range(mod.getNumParameters()):
    element = mod.getParameter(i)
    generated_code.extend('  {0} = {1}\n'.format(element.getId(), element.getValue()))
  
  # write out boundary species 
  generated_code.extend('  \n')
  generated_code.extend('  #boundary species\n')
  for i in range(mod.getNumSpecies()):
    element = mod.getSpecies(i)
    if element.getBoundaryCondition() == False: 
      continue 
    if element.isSetInitialConcentration(): 
      generated_code.extend('  {0} = {1}\n'.format(element.getId(), element.getInitialConcentration()))
    else:
      generated_code.extend('  {0} = {1}\n'.format(element.getId(), element.getInitialAmount()))  
  
  # write derive function
  generated_code.extend('  \n')
  generated_code.extend('  def ode_fun(__Y__, t):\n')
  for i in range(len(variables.keys())): 
    generated_code.extend('    {0} = __Y__[{1}]\n'.format(variables.keys()[i], i))
  generated_code.extend('\n')
  
  for i in range(mod.getNumReactions()): 
    reaction = mod.getReaction(i)
    kinetics = reaction.getKineticLaw()  
    
    generated_code.extend('    {0} = {1}\n'.format(reaction.getId(),  kinetics.getFormula()))
    
    for j in range(reaction.getNumReactants()): 
      ref = reaction.getReactant(j)
      species = mod.getSpecies(ref.getSpecies())
      if species.getBoundaryCondition() == True: 
        continue
      if ref.getStoichiometry() == 1.0: 
        variables[species.getId()].append('-{0}'.format(reaction.getId()))
      else:
        variables[species.getId()].append('-({0})*{1}'.format(ref.getStoichiometry(), reaction.getId()))
    for j in range(reaction.getNumProducts()): 
      ref = reaction.getProduct(j)
      species = mod.getSpecies(ref.getSpecies())
      if species.getBoundaryCondition() == True: 
        continue
      if ref.getStoichiometry() == 1.0: 
        variables[species.getId()].append('{0}'.format(reaction.getId()))
      else:
        variables[species.getId()].append('({0})*{1}'.format(ref.getStoichiometry(), reaction.getId()))
  
  generated_code.extend('\n')
    
  generated_code.extend('    return array([')
  for i in range(len(variables.keys())):
    for eqn in variables[variables.keys()[i]]:
      generated_code.extend(' + ({0})'.format(eqn))
    if i + 1 < len(variables.keys()):
      generated_code.extend(',\n      ')
  generated_code.extend('    ])\n')
  generated_code.extend('\n')
  
  generated_code.extend('  time = linspace(t0, tend, numPoints)\n')
  
  # 
  # write out initial concentrations 
  # 
  generated_code.extend('  yinit= array([')
  count = 0
  for key in variables.keys():
    # get initialValue 
    element = mod.getElementBySId(key)
    if element.getTypeCode() == SBML_PARAMETER: 
      generated_code.extend('{0}'.format(element.getValue()))
    elif element.getTypeCode() == SBML_SPECIES: 
      if element.isSetInitialConcentration(): 
        generated_code.extend('{0}'.format(element.getInitialConcentration()))
      else: 
        generated_code.extend('{0}'.format(element.getInitialAmount()))
    else: 
      generated_code.extend('{0}'.format(element.getSize()))
    count += 1
    if count < len(variables.keys()):
      generated_code.extend(', ')
  generated_code.extend('])\n')
  generated_code.extend('  \n')
  
  generated_code.extend('  y = odeint(ode_fun, yinit, time)\n')
  generated_code.extend('\n')
  
  generated_code.extend('  return time, y\n')
  generated_code.extend('\n')
  generated_code.extend('\n')
  generated_code.extend('time, result = simulateModel({0}, {1}, {2})\n'.format(t0, tEnd, numPoints))
  generated_code.extend('\n')
  
  #
  # write out plotting code 
  #
  generated_code.extend('fig = figure()\n')
  generated_code.extend('ax = subplot(111)\n')
  
  for i in range(len(variables.keys())): 
    generated_code.extend('plot(time,result[:,{0}], label="{1}", lw=1.5)\n'.format(i, variables.keys()[i]))
    
  generated_code.extend('box = ax.get_position()\n')
  generated_code.extend('ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])\n')
  generated_code.extend('xlabel("time")\n')
  generated_code.extend('ylabel("concentration")\n')
  generated_code.extend('legend(loc="center left", bbox_to_anchor=(1, 0.5))\n')
  generated_code.extend('show()\n')
  
  
  # convert generated code to string 
  result = str(generated_code);
  return result 
  
  
if __name__ == "__main__":
  result = generateCodeForFile(sys.argv[1])

  file = open('generated.py', 'w')
  file.write(result)
  file.close()
  
  exec(result)

  