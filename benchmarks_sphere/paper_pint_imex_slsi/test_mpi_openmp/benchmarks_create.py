#! /usr/bin/env python3
#
#  Create series of job to be run with sweet
#
#  Pedro Peixoto <pedrosp@ime.usp.br>
#  modified from Martin Schreiber initial job_create
#
#-------------------------------------------------------

import os
import sys
import stat
import math
from glob import glob

#Classes containing sweet compile/run basic option
from mule_local.JobGeneration import *
from mule_local.SWEETRuntimeParametersScenarios import *
from mule.JobParallelization import *
from mule.JobParallelizationDimOptions import *

sim_type = sys.argv[1];
if sim_type == "xbraid":
    folder_fine_sim = "fine_sim";

    ## find fine, reference simulation
    a = [name for name in os.listdir(folder_fine_sim)
            if os.path.isdir(os.path.join(folder_fine_sim, name))]
    assert len(a) == 1;
    ref_job = os.path.abspath(os.getcwd()) + "/" + folder_fine_sim + "/" + a[0];


#Create main compile/run options
jg = JobGeneration()

#Get Earth parameters (if necessary)
earth = EarthMKSDimensions()

#
# Run simulation on plane or sphere
#
#Basic plane options
jg.compile.program = "swe_sphere"
jg.compile.mode = "debug"
jg.compile.sweet_mpi = "enable"


jg.compile.sphere_spectral_space = 'enable';
jg.compile.sphere_spectral_dealiasing = 'enable';

# Verbosity mode
jg.runtime.verbosity = 3

jg.runtime.output_file_mode = 'bin';

#
# Benchmark ID
# 14: Steady diagonal benchmark
#
#jg.runtime.bench_id = 1
jg.runtime.benchmark_name = "three_gaussian_bumps_phi_pint"

#
# Compute error or difference to initial data
#
####jg.runtime.compute_error = 0

# Enable/Disbale GUI
jg = DisableGUI(jg)

#
# REXI method
jg.runtime.rexi_method = 'direct'
#jg.runtime.rexi_use_direct_solution = 1

# Parameters for SL-REXI paper
#-----------------------------
jg = RuntimeSWEPlaneEarthParam(jg)
#jg = RuntimeSWENonDimParam(jg)

jg.runtime.viscosity = 0


###max_simulation_time = 102400.;
max_simulation_time = 60 * 60 * 2;

#
# Time, Mode and Physical resolution
#
timestep_size_reference = 60.;
timestep_size_fine = 60.; #3600 #1 hour  #864000/10 #1 day

jg.runtime.max_simulation_time = max_simulation_time; #1 day #timestep_size_reference #864000 #10 days
jg.runtime.output_timestep_size = max_simulation_time;

jg.runtime.timestep_size = timestep_size_reference
jg.runtime.timestepping_method = "l_irk_n_erk"
jg.runtime.timestepping_order = 2
jg.runtime.timestepping_order2 = 2
jg.runtime.space_res_physical = -1
jg.runtime.space_res_spectral = 256
jg.runtime.exp_direct_precompute_phin = 1;

if sim_type == "xbraid":
    jg.compile.xbraid = "mpi";
    jg.runtime.xbraid_enabled = 1;
    jg.runtime.xbraid_min_coarse = 2
    jg.runtime.xbraid_nrelax0 = -1
    jg.runtime.xbraid_tol = 0.
    jg.runtime.xbraid_tnorm = 2
    jg.runtime.xbraid_cfactor0 = -1
    jg.runtime.xbraid_max_iter = 11
    jg.runtime.xbraid_res = 0
    jg.runtime.xbraid_storage = 0
    jg.runtime.xbraid_print_level = 2
    jg.runtime.xbraid_access_level = 2
    jg.runtime.xbraid_run_wrapper_tests = 0
    jg.runtime.xbraid_fullrnorm = 2
    jg.runtime.xbraid_use_seq_soln = 0
    jg.runtime.xbraid_use_rand = 1
    jg.runtime.xbraid_pt = 1
    jg.runtime.xbraid_timestepping_method = "l_irk_n_erk,lg_exp_na_sl_lc_nr_etdrk_uv"
    jg.runtime.xbraid_timestepping_order = 2
    jg.runtime.xbraid_timestepping_order2 = 2
    jg.runtime.xbraid_verbosity = 0;
    jg.runtime.xbraid_viscosity_order = 2;
    jg.runtime.xbraid_viscosity_coefficient = str(jg.runtime.viscosity) + ",1e6";

    jg.runtime.xbraid_load_ref_csv_files = 1;
    jg.runtime.xbraid_path_ref_csv_files = ref_job;
    jg.runtime.xbraid_load_fine_csv_files = 1;
    jg.runtime.xbraid_path_fine_csv_files = ref_job;
    jg.runtime.xbraid_store_iterations = 0;
    #######jg.runtime.xbraid_spectral_ref = 1024;

    jg.runtime.xbraid_timestepping_method = "l_irk_n_erk,l_irk_n_erk";
    jg.runtime.xbraid_cfactor = 2;
    jg.runtime.xbraid_max_levels = 2;
    jg.runtime.xbraid_skip = 1;
    jg.runtime.timestep_size = timestep_size_reference;
    jg.runtime.xbraid_spatial_coarsening = 51;
    jg.runtime.xbraid_fmg = 1;
    jg.runtime.xbraid_fmg_vcyc = 1;
    jg.runtime.xbraid_nrelax = 0;
    
    nb_pts = [1, 2, 4]; ## number of parallel processors in time


    for nb_pt in nb_pts:

        jg.runtime.xbraid_pt = nb_pt;

        params_pspace_num_cores_per_rank = [jg.platform_resources.num_cores_per_socket]
        #params_pspace_num_threads_per_rank = [i for i in range(1, jg.platform_resources.num_cores_per_socket+1)]
        params_pspace_num_threads_per_rank = [jg.platform_resources.num_cores_per_socket]
        params_ptime_num_cores_per_rank = [1]

        # Update TIME parallelization
        ptime = JobParallelizationDimOptions('time')
        ptime.num_cores_per_rank = 1
        ptime.num_threads_per_rank = 1 #pspace.num_cores_per_rank
        ptime.num_ranks = nb_pt

        pspace = JobParallelizationDimOptions('space')
        pspace.num_cores_per_rank = 32
        ###pspace.num_threads_per_rank = params_pspace_num_cores_per_rank[-1]
        pspace.num_threads_per_rank = 16
        pspace.num_ranks = 1

        # Setup parallelization
        jg.setup_parallelization([pspace, ptime])


        jg.gen_jobscript_directory()

elif sim_type == "ref":
####### fine simulation
    jg.compile.parareal = "none";
    jg.compile.xbraid = "none";
    jg.runtime.parareal_enabled = 0;
    jg.runtime.xbraid_enabled = 0;


    params_pspace_num_cores_per_rank = [jg.platform_resources.num_cores_per_socket]
    #params_pspace_num_threads_per_rank = [i for i in range(1, jg.platform_resources.num_cores_per_socket+1)]
    params_pspace_num_threads_per_rank = [jg.platform_resources.num_cores_per_socket]
    params_ptime_num_cores_per_rank = [1]

    # Update TIME parallelization
    ptime = JobParallelizationDimOptions('time')
    ptime.num_cores_per_rank = 1
    ptime.num_threads_per_rank = 1 #pspace.num_cores_per_rank
    ptime.num_ranks = 1

    pspace = JobParallelizationDimOptions('space')
    pspace.num_cores_per_rank = 32
    ###pspace.num_threads_per_rank = params_pspace_num_cores_per_rank[-1]
    pspace.num_threads_per_rank = 16
    pspace.num_ranks = 1

    # Setup parallelization
    jg.setup_parallelization([pspace, ptime])

    jg.parallelization.mpiexec_disabled = False
    ####jg.parallelization.mpiexec_disabled = True

    jg.gen_jobscript_directory();
