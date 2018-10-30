#PySWMM test code

import pyswmm as ps                              # import package

sim = ps.Simulation('tutorial1.inp')             # load & run simulation with no interaction
sim.execute()

with ps.Simulation('tutorial1.inp') as sim:      # load and run simulation and print each timestep
	for step in sim:
		print(sim.current_time)
	sim.report()
	