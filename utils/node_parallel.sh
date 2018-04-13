#!/bin/bash
# A cross-node version of gnu parallel
#
#    node-parallel [command [arguments]] < list_of_arguments
#    node-parallel [command [arguments]] ( ::: arguments | :::: argfile(s) )*
#
# The command can use gnu parallel replace tags (e.g. '{}' is replaced with an
# argument), see 'man parallel' for full details.
# Note that unlike gnu parallel a '{}' will not be automatically added if not
# present, it must be manually added.
#
# The command is executed in $PBS_NCPUS parallel jobs, on all the cpus
# available
#
# Set the environment variable $NP_DEBUG to a non-empty string to show the
# commands being executed

set -eu

# Save the environment in the shared tmp directory so we can see it from all nodes
NP_ENVFILE=$(mktemp --tmpdir="/short/$PROJECT/$USER/tmp")
NP_ERRFILE=$(mktemp --tmpdir="/short/$PROJECT/$USER/tmp")
env | sed -e '/{/,/}/d' -e 's/^.*$/export '"'"'\0'"'"'/' > $NP_ENVFILE

# Split arguments into command and inputs (empty argument at the end  handles
# the case where there are no inputs)
ARGS="$* :::"
COMMAND="${ARGS%%:::*}"
INPUTS="${ARGS#*:::}"
INPUTS="${INPUTS:+:::}${INPUTS}"

# Load required modules
module load parallel
module load pbs

# We run ${PBS_NCPUS} jobs in parallel
# ... Each on its own CPU
# ... Loading the enviornment from the parent job

PBSDSH=${PBS_NCPUS:+"pbsdsh -n {%} --"}

parallel ${NP_DEBUG:+"--verbose"} --jobs ${PBS_NCPUS:-2} ${NP_PARALLEL_ARGS} -- \
    ${PBSDSH} \
    /bin/bash -l -c "'\
        source ${NP_ENVFILE}; \
        (${COMMAND}); \
        echo \$? >> ${NP_ERRFILE}'" \
    ${INPUTS}

# pbsdsh doesn't propogate errors, so collect them manually
ERR=$(grep -v '0' ${NP_ERRFILE} | wc -l)

# Remove temporary files
rm ${NP_ENVFILE}
rm ${NP_ERRFILE}

exit ${ERR}
