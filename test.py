##########################################################################
# USEAGE: python test.py forcefield_type True/False xml_file
# forcefield_type can be tip5p or mbpol
# xml_file can be ignored and will default to tip5p.xml or mbpol.xml 
#    depending on the forcefield_type
##########################################################################

from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys

forcefield_type = sys.argv[1]
if sys.argv[2] == 'True':
    rigid_water = True
else:
    rigid_water = False


if forcefield_type == 'tip5p':
    water = 'water3_tip5p.pdb'
else:
    water = 'water3.pdb'



pdb = app.PDBFile(water)
#pdb = app.PDBFile('a.pdb')
#pdb = app.PDBFile('water14_cluster.pdb')
if len(sys.argv) == 4:
    xml_file = sys.argv[3]
elif forcefield_type == 'tip5p':
    xml_file = 'tip5p.xml'
else:
    xml_file = '/home/sebastian/anaconda3/lib/python3.4/site-packages/mbpol-1.0-py3.4-linux-x86_64.egg/mbpol.xml'

forcefield = app.ForceField(xml_file)

print ('Running with {} {}'.format(water, xml_file))
    
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic, 
    nonbondedCutoff=0.9*unit.nanometers, constraints=None, rigidWater=rigid_water,
    ewaldErrorTolerance=0.0005)
integrator = mm.VerletIntegrator(0.2*unit.femtoseconds)

platform = mm.Platform.getPlatformByName('Reference')
#platform = mm.Platform.getPlatformByName('CUDA')
#properties = {'CudaPrecision': 'mixed'}
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.context.computeVirtualSites()


simulation.context.setVelocitiesToTemperature(300*unit.kelvin)

# simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True, time=True, 
#    potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, 
#    progress=True, remainingTime=True, speed=True, totalSteps=1100, 
#    separator='\t'))

filename = forcefield_type + '_trajectory.pdb'
if rigid_water == False:
    filename = 'flexible_' + filename

print ('Writing to {}'.format(filename))

simulation.reporters.append(app.PDBReporter(filename, 1))

print('Running Production...')
simulation.step(100)
print('Done!')

