
/*
 * main.cpp
 *
 * PFASST SWE on the sphere implementation
 *
 *  Created on: 30 Nov 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#include <sweet/FatalError.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/SphereDiagnostics.hpp>
#include <benchmarks_sphere/SphereBenchmarksCombined.hpp>

#include "libpfasst_swe_sphere/LevelSingleton.hpp"
#include "libpfasst_swe_sphere/SphereDataCtx.hpp"
#include <mpi.h>

#define WITH_MPI

SimulationVariables simVars;
std::vector<LevelSingleton> levelSingletons;

extern "C"
{
   /* Driver function for pfasst control */
  void fmain (SphereDataCtx* pd_ctx, 
	      const int*     nlevels, 
	      const int*     niters, 
	      const int*     nsweeps_coarse,
	      const int      nnodes[], 
	      const char*    qtype_name, 
	      const int*     qtype_name_len,
 	      const int*     nfields, 
	      const int      nvars_per_field[], 
 	      double*        t_max, 
	      double*        dt);
}

/**
 * Main function launching LibPFASST
 */

int main(int i_argc, char *i_argv[])
{
  MPI_Init(&i_argc, &i_argv);

  // input parameter names (specific ones for this program)
  const char *bogus_var_names[] = {
    "rexi-use-coriolis-formulation",
    "compute-error",
    nullptr
  };

  // default values for specific input (for general input see SimulationVariables.hpp)
  simVars.bogus.var[0] = 1;
  simVars.bogus.var[1] = 1;
  
  // Help menu
  if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
    {
      std::cout << "--compute-error [0/1]Output errors (if available, default: 1)" << std::endl;
      std::cout << "--rexi-use-coriolis-formulation [0/1]Use Coriolisincluding  solver for REXI (default: 1)" << std::endl;
      return -1;
    }

  simVars.timecontrol.current_timestep_size = - simVars.sim.CFL; 

  // define the number of levels and SDC nodes for each level
  // note: level #nlevels-1 is the finest, level #0 is the coarsest

  int nnodes[simVars.libpfasst.nlevels];
  nnodes[simVars.libpfasst.nlevels-1] = simVars.libpfasst.nnodes; // finest level

  switch (simVars.libpfasst.nlevels) 
    {
      // One level (nothing to do)
      case 1: {
	break;
      }
      // Two levels
      case 2: {
	if (simVars.libpfasst.nnodes == 3)
	  nnodes[0] = 2;
	else if (simVars.libpfasst.nnodes == 5 || 
		 simVars.libpfasst.nnodes == 9)
	  nnodes[0] = 3; 
	else 
	  FatalError("With 2 levels, the number of SDC nodes on the fine level must be either 3, 5, or 9");
	break;
      }
      // Three levels
      case 3: {
	if (simVars.libpfasst.nnodes == 9) 
	  {
	    nnodes[0] = 3; 
	    nnodes[1] = 5;
	  }
	else if (simVars.libpfasst.nnodes == 5) 
	  {
	    nnodes[0] = 2; 
	    nnodes[1] = 3;
	  }
	else if (simVars.libpfasst.nnodes == 3) 
	  {
	    nnodes[0] = 3; 
	    nnodes[1] = 3;
	  }
	else 
	  FatalError("With 3 levels, the number of SDC nodes on the fine level must be either 5, or 9");
	break;
      }
      // All other cases not supported yet
      default:
	FatalError("Only 1, 2, or 3 levels are currently supported");
    }
  
  // setup the LevelSingletons for all levels
  // note: level #nlevels-1 is the finest, level #0 is the coarsest

  levelSingletons.resize(simVars.libpfasst.nlevels);

  // setup the finest level singleton
  const int fineLevelId = simVars.libpfasst.nlevels-1;


  levelSingletons[fineLevelId].level = fineLevelId;

  // setup data configuration in fine level

  levelSingletons[fineLevelId].dataConfig.setupAuto(
						    simVars.disc.res_physical,
						    simVars.disc.res_spectral 
						    );
  
  int res_physical_nodealiasing[2] = {
    2*(simVars.disc.res_spectral[0]+1),
    simVars.disc.res_spectral[1]+2
  };

  levelSingletons[fineLevelId].dataConfigNoDealiasing.setupAuto(
								res_physical_nodealiasing,
								simVars.disc.res_spectral 
								);

  // setup data operators in fine level
  
  levelSingletons[fineLevelId].op.setup(
					&(levelSingletons[fineLevelId].dataConfig),
					simVars.sim.earth_radius
					);
  levelSingletons[fineLevelId].opNoDealiasing.setup(
						    &(levelSingletons[fineLevelId].dataConfigNoDealiasing),
						    simVars.sim.earth_radius
						    );
  
  // define the number of modes for the coarser levels
  for (int i = 1; i < simVars.libpfasst.nlevels; i++)
    {
      const int thisLevelId = simVars.libpfasst.nlevels-1-i;
      levelSingletons[thisLevelId].level = thisLevelId;
	    
      // setup data configuration at this level

      levelSingletons[thisLevelId].dataConfig.setupAdditionalModes(
								   &(levelSingletons[simVars.libpfasst.nlevels-i].dataConfig),
								   -std::ceil(simVars.disc.res_spectral[0]*pow(simVars.libpfasst.coarsening_multiplier,i)),
								   -std::ceil(simVars.disc.res_spectral[1]*pow(simVars.libpfasst.coarsening_multiplier,i))
								   );
      
      // setup data operators at this level
  
      levelSingletons[thisLevelId].op.setup(
					    &(levelSingletons[thisLevelId].dataConfig),
					    simVars.sim.earth_radius
					    );
    }

  // define the SWEET parameters

  const int nfields = 3;  // number of vector fields (here, height and two horizontal velocities)
  int nvars_per_field[simVars.libpfasst.nlevels];
  for (int i = 0; i < simVars.libpfasst.nlevels; ++i) 
    nvars_per_field[i] = 2*levelSingletons[i].dataConfig.spectral_array_data_number_of_elements;  // number of degrees of freedom per vector field

  // initialize the topography before instantiating the SphereDataCtx object
  if (simVars.setup.benchmark_scenario_name == "flow_over_mountain") 
  {

      // create h_topo with the configuration at the finest level
      simVars.sim.h_topo.setup(&(levelSingletons[simVars.libpfasst.nlevels-1].dataConfig));

      // initialize the topography
      SphereBenchmarksCombined::setupTopography(simVars,
						levelSingletons[simVars.libpfasst.nlevels-1].op
						);
  }

  // instantiate the SphereDataCtx object 
  SphereDataCtx* pd_ctx = new SphereDataCtx(
					    &simVars,
					    &levelSingletons,
					    nnodes
					    );

  // get the C string length (needed by Fortran...)
  int string_length = simVars.libpfasst.nodes_type.size();

  // call LibPFASST to advance in time
  fmain(
   	pd_ctx,                                       // user defined context
	&simVars.libpfasst.nlevels,                   // number of SDC levels
	&simVars.libpfasst.niters,                    // number of SDC iterations
	&simVars.libpfasst.nsweeps_coarse,            // number of SDC sweeps on coarse level
	nnodes,                                       // number of SDC nodes 
 	(simVars.libpfasst.nodes_type).c_str(),       // type of nodes
	&string_length,                               // length of (simVars.libpfasst.nodes_type).c_str()
	&nfields,                                     // number of vector fields
	nvars_per_field,                              // number of dofs per vector field
 	&(simVars.timecontrol.max_simulation_time),   // simulation time
   	&(simVars.timecontrol.current_timestep_size)  // time step size
 	); 

  // release the memory
  delete pd_ctx;

  MPI_Finalize();
}

