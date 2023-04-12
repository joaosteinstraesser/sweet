#! /bin/bash

cd "$(dirname $0)"


TIMESTEPPING_GROUP="l1"
REXI_PHI_PRECOMP="1"

COMMON="../swe_plane_timestepper_convergence_common_no_test/"


mule.benchmark.cleanup_all || exit 1

$COMMON/benchmark_create_job_scripts.py $TIMESTEPPING_GROUP $REXI_PHI_PRECOMP || exit 1

mule.benchmark.jobs_run_directly || exit 1

$COMMON/postprocessing_pickle.py || exit 1

$COMMON/postprocessing_convergence_test.py || exit 1

mule.benchmark.cleanup_all || exit 1
