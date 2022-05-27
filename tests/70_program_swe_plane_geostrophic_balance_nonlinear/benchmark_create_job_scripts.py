#! /usr/bin/env python3

import sys
from itertools import product

from mule_local.JobMule import *
from mule.exec_program import *
from mule.InfoError import *
from mule_local.SWEETRuntimeParametersScenarios import *

jg = JobGeneration()


"""
Compile parameters
"""
###params_compile_sweet_mpi = ['enable', 'disable']
###params_compile_threading = ['omp', 'off']
###params_compile_thread_parallel_sum = ['enable', 'disable']
params_compile_sweet_mpi = ['disable']
params_compile_threading = ['omp']
params_compile_thread_parallel_sum = ['disable']

jg = CompileSWEPlane(jg)
jg = RuntimeSWEPlaneEarthParam(jg)

jg.compile.plane_spectral_space = 'enable'
jg.compile.plane_spectral_dealiasing = 'enable'

jg.compile.quadmath = 'enable'
jg.unique_id_filter = ['runtime.simparams', 'parallelization', 'benchmark', 'runtime.rexi_params']


"""
Runtime parameters
"""
##params_runtime_timestep_sizes = [30]
params_runtime_timestep_sizes = [3600]

jg.runtime.benchmark_name = 'steady_state_zonal_flow'
##jg.runtime.benchmark_name = 'rotated_steady_state'

jg.runtime.space_res_spectral = 128
jg.runtime.space_res_physical = None

jg.runtime.compute_error = 1

# run 10 time steps
##jg.runtime.max_simulation_time = 30*10
jg.runtime.max_simulation_time = 864000


# Use moderate CI-REXI values
# Set later on
#jg.runtime.rexi_method = 'ci'
#jg.runtime.rexi_ci_n = 16
#jg.runtime.rexi_ci_max_real = 1
#jg.runtime.rexi_ci_max_imag = 1
#jg.runtime.rexi_ci_mu = 0
#jg.runtime.rexi_ci_primitive = 'circle'
#jg.runtime.rexi_sphere_preallocation = 1
jg.runtime.rexi_method = 'direct'

jg.runtime.instability_checks = 0
jg.runtime.verbosity = 10

# output results after end
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time
jg.runtime.output_filename = "-"




"""
Parallelization parameters
"""

# Update TIME parallelization
ptime = JobParallelizationDimOptions('time')
ptime.num_cores_per_rank = 1
ptime.num_threads_per_rank = 1

if jg.platform_resources.num_cores_per_node <= 1:
    ptime.num_ranks = 1
else:
    ptime.num_ranks = 2


pspace = JobParallelizationDimOptions('space')
pspace.num_cores_per_rank = 1
pspace.num_threads_per_rank = jg.platform_resources.num_cores_per_socket//ptime.num_threads_per_rank
pspace.num_ranks = 1

if pspace.num_threads_per_rank == 0:
    pspace.num_threads_per_rank = 1

# Setup parallelization
jg.setup_parallelization([pspace, ptime])




ts_methods = [
    ['ln_erk',        4,    4,    0],    # reference solution

    ##['ln_erk',        2,    2,    0],

    ##['l_erk_n_erk',        2,    2,    0],

    ['l_irk_n_erk',        1,    1,    0],

    ['l_cn_na_sl_nd_settls',        2,    2,    0],

    ['l_rexi_na_sl_nd_etdrk',        1,    1,    0],
    ['l_rexi_na_sl_nd_etdrk',        2,    2,    0],
    ['l_rexi_na_sl_nd_etdrk',        4,    4,    0],
    ['l_rexi_na_sl_nd_etdrk',        -2,    -2,    0],

    ['l_rexi_na_sl_nd_settls',    2,    2,    0],

    ['l_rexi_n_etdrk',    2,    2,    0],



]


interp_orders = [1, 5, 32, 34, 36, 38, -3]
coriolis_treatments = ['linear', 'nonlinear'];

#
# allow including this file
#
if __name__ == "__main__":

    #
    # Create job scripts
    #
    for tsm in ts_methods[1:]:

        for interp_order in interp_orders:

            for coriolis_treatment in coriolis_treatments:

                if (not "sl" in tsm[0]) and not (interp_order == 1):
                    continue;

                if (not "sl" in tsm[0]) and not (coriolis_treatment == 'linear'):
                    continue;

                jg.runtime.timestepping_method = tsm[0]
                jg.runtime.timestepping_order = tsm[1]
                jg.runtime.timestepping_order2 = tsm[2]
                jg.runtime.semi_lagrangian_interpolation_order = interp_order;
                jg.runtime.coriolis_treatment = coriolis_treatment;

                if len(tsm) > 4:
                    s = tsm[4]
                    jg.runtime.load_from_dict(tsm[4])

                for jg.runtime.timestep_size in params_runtime_timestep_sizes:

                    for (
                        jg.compile.threading,
                        jg.compile.rexi_thread_parallel_sum,
                        jg.compile.sweet_mpi
                    ) in product(
                        params_compile_threading,
                        params_compile_thread_parallel_sum,
                        params_compile_sweet_mpi
                    ):
                        if 'exp_' in jg.runtime.timestepping_method:

                            jg.runtime.rexi_method = 'direct'
                            jg.gen_jobscript_directory()
                            jg.runtime.rexi_method = ''

                        else:
                            if jg.compile.sweet_mpi == 'enable':
                                    continue

                            if jg.compile.rexi_thread_parallel_sum == 'enable':
                                    continue

                            jg.gen_jobscript_directory()


