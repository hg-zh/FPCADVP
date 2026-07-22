#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=570G
#SBATCH --time=120:00:00
#SBATCH --job-name=sensitivity

set -Eeuo pipefail
umask 027

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "ERROR: this script must run inside a Slurm allocation." >&2
  exit 2
fi

code_dir_input="${SENSITIVITY_CODE_DIR:-${SLURM_SUBMIT_DIR:-}}"
if [[ -z "${code_dir_input}" || ! -d "${code_dir_input}" ]]; then
  echo "ERROR: set SENSITIVITY_CODE_DIR to the sensitivity_final directory." >&2
  exit 3
fi
readonly CODE_DIR="$(cd "${code_dir_input}" && pwd -P)"

output_root_input="${SENSITIVITY_OUTPUT_ROOT:-${CODE_DIR}}"
log_root_input="${SENSITIVITY_LOG_ROOT:-${CODE_DIR}/logs}"
mkdir -p "${output_root_input}" "${log_root_input}"
readonly OUTPUT_ROOT="$(cd "${output_root_input}" && pwd -P)"
readonly RESULTS_ROOT="${OUTPUT_ROOT}/results"
readonly LOG_ROOT="$(cd "${log_root_input}" && pwd -P)"
mkdir -p "${RESULTS_ROOT}"

if [[ ! -w "${RESULTS_ROOT}" || ! -w "${LOG_ROOT}" ]]; then
  echo "ERROR: results and log directories must be writable." >&2
  exit 3
fi

module_input="${TRANSPORT_FPCA_MODULE:-${CODE_DIR}/../Simulation/parallel_run/transport_fpca_modular.R}"
if [[ "${module_input}" != /* ]]; then
  module_input="${CODE_DIR}/${module_input}"
fi
if [[ ! -r "${module_input}" ]]; then
  echo "ERROR: transport-FPCA dependency is missing or unreadable: ${module_input}" >&2
  exit 3
fi
readonly MODULE_SCRIPT="$(cd "$(dirname "${module_input}")" && pwd -P)/$(basename "${module_input}")"

for required_file in run_sensitivity.R aggregate_results.R generate_table.R run_sensitivity.sbatch; do
  if [[ ! -r "${CODE_DIR}/${required_file}" ]]; then
    echo "ERROR: required runtime file is missing or unreadable: ${CODE_DIR}/${required_file}" >&2
    exit 3
  fi
done

readonly RUN_TYPE="${RUN_TYPE:-full}"
case "${RUN_TYPE}" in
  smoke)
    readonly REPS=2
    readonly N_VEC="20"
    readonly OBS_VEC="3"
    readonly M_VEC="200"
    readonly EXPECTED_ROWS=30
    readonly WRITE_TABLE=FALSE
    ;;
  full)
    readonly REPS=200
    readonly N_VEC="20,50,100,200,400"
    readonly OBS_VEC="3,5,10,20"
    readonly M_VEC="200"
    readonly EXPECTED_ROWS=60000
    readonly WRITE_TABLE=TRUE
    ;;
  *)
    echo "ERROR: RUN_TYPE must be 'smoke' or 'full'; got '${RUN_TYPE}'." >&2
    exit 2
    ;;
esac

readonly RUN_ID="${RUN_TYPE}_${SLURM_JOB_ID}"
readonly RUN_DIR="${RESULTS_ROOT}/${RUN_ID}"
if [[ -e "${RUN_DIR}" ]]; then
  echo "ERROR: refusing to overwrite existing run directory: ${RUN_DIR}" >&2
  exit 3
fi

RUN_DIR_OWNED=0
COMPLETION_KIND="SUCCESS"
finish() {
  local status=$?
  trap - EXIT
  if (( RUN_DIR_OWNED == 1 )); then
    mkdir -p "${RUN_DIR}"
    if (( status == 0 )); then
      if [[ "${COMPLETION_KIND}" == "FIT_FAILURES" ]]; then
        printf 'run_id=%s\njob_id=%s\nrun_type=%s\nstatus=complete_with_fit_failures\ncompleted_utc=%s\n' \
          "${RUN_ID}" "${SLURM_JOB_ID}" "${RUN_TYPE}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
          > "${RUN_DIR}/_COMPLETE_WITH_FIT_FAILURES"
      else
        printf 'run_id=%s\njob_id=%s\nrun_type=%s\nstatus=complete\ncompleted_utc=%s\n' \
          "${RUN_ID}" "${SLURM_JOB_ID}" "${RUN_TYPE}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
          > "${RUN_DIR}/_SUCCESS"
      fi
    else
      printf 'run_id=%s\njob_id=%s\nrun_type=%s\nexit_status=%s\nfailed_utc=%s\n' \
        "${RUN_ID}" "${SLURM_JOB_ID}" "${RUN_TYPE}" "${status}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
        > "${RUN_DIR}/_FAILED"
    fi
  fi
  exit "${status}"
}
trap finish EXIT

RUN_DIR_OWNED=1
mkdir -p "${RUN_DIR}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

if [[ -n "${R_MODULE:-}" ]]; then
  if ! type module >/dev/null 2>&1; then
    echo "ERROR: R_MODULE is set but the environment-modules command is unavailable." >&2
    exit 4
  fi
  module load "${R_MODULE}"
fi
if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript is unavailable." >&2
  exit 4
fi
if ! command -v sha256sum >/dev/null 2>&1; then
  echo "ERROR: sha256sum is required for the run manifest." >&2
  exit 4
fi
if [[ ! "${SLURM_CPUS_PER_TASK:-}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: SLURM_CPUS_PER_TASK is missing or invalid." >&2
  exit 4
fi

Rscript -e '
  expected_r <- "4.5.0"
  expected_fdapace <- "0.6.0"
  actual_r <- as.character(getRversion())
  if (!identical(actual_r, expected_r)) {
    stop("R version mismatch: expected ", expected_r, ", got ", actual_r)
  }
  if (!requireNamespace("fdapace", quietly = TRUE)) {
    stop("R package fdapace is unavailable")
  }
  actual_fdapace <- as.character(utils::packageVersion("fdapace"))
  if (!identical(actual_fdapace, expected_fdapace)) {
    stop("fdapace version mismatch: expected ", expected_fdapace,
         ", got ", actual_fdapace)
  }
'

cd "${CODE_DIR}"
Rscript run_sensitivity.R \
  mode="${RUN_TYPE}" \
  reps="${REPS}" \
  n_vec="${N_VEC}" \
  N_vec="${OBS_VEC}" \
  m_vec="${M_VEC}" \
  kernels=gauss,epan,quar \
  bw_multipliers=0.5,0.75,1,1.25,1.5 \
  cores="${SLURM_CPUS_PER_TASK}" \
  seed_base=2000 \
  setting=1 \
  decayR=1 \
  latent=FALSE \
  methodXi=CE \
  t_design=random \
  t_grid_len=51 \
  eps_u=0.05 \
  u_cap=0.95 \
  Z_cap=8 \
  beta_a=0.6,3 \
  beta_b=0.6,3 \
  mod_script="${MODULE_SCRIPT}" \
  check_only=TRUE

readonly BATCH_STARTED_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
readonly DRIVER_SHA256="$(sha256sum "${CODE_DIR}/run_sensitivity.R" | awk '{print $1}')"
readonly AGGREGATOR_SHA256="$(sha256sum "${CODE_DIR}/aggregate_results.R" | awk '{print $1}')"
readonly TABLE_SHA256="$(sha256sum "${CODE_DIR}/generate_table.R" | awk '{print $1}')"
readonly BATCH_SHA256="$(sha256sum "${CODE_DIR}/run_sensitivity.sbatch" | awk '{print $1}')"
readonly MODULE_SHA256="$(sha256sum "${MODULE_SCRIPT}" | awk '{print $1}')"
ACTUAL_TIME_LIMIT="$(squeue -h -j "${SLURM_JOB_ID}" -o '%l' 2>/dev/null || true)"
if [[ -z "${ACTUAL_TIME_LIMIT}" ]]; then
  ACTUAL_TIME_LIMIT="${SLURM_TIMELIMIT:-<unavailable>}"
fi
readonly ACTUAL_TIME_LIMIT

cat > "${RUN_DIR}/batch_manifest.txt" <<EOF
run_id=${RUN_ID}
run_type=${RUN_TYPE}
job_id=${SLURM_JOB_ID}
job_name=${SLURM_JOB_NAME:-sensitivity}
partition=${SLURM_JOB_PARTITION:-<unspecified>}
account=${SLURM_JOB_ACCOUNT:-<unspecified>}
submit_dir=${SLURM_SUBMIT_DIR:-<unset>}
node_list=${SLURM_JOB_NODELIST:-<unset>}
nodes=${SLURM_NNODES:-1}
tasks=${SLURM_NTASKS:-1}
cpus_per_task=${SLURM_CPUS_PER_TASK}
memory_per_node_mb=${SLURM_MEM_PER_NODE:-<unset>}
time_limit=${ACTUAL_TIME_LIMIT}
log_root=${LOG_ROOT}
rscript=$(command -v Rscript)
started_utc=${BATCH_STARTED_UTC}
expected_rows=${EXPECTED_ROWS}
reps=${REPS}
n_vec=${N_VEC}
N_vec=${OBS_VEC}
m_vec=${M_VEC}
kernels=gauss,epan,quar
bw_multipliers=0.5,0.75,1,1.25,1.5
run_sensitivity_sha256=${DRIVER_SHA256}
aggregate_results_sha256=${AGGREGATOR_SHA256}
generate_table_sha256=${TABLE_SHA256}
run_sensitivity_sbatch_sha256=${BATCH_SHA256}
transport_fpca_modular_sha256=${MODULE_SHA256}
EOF

echo "Starting ${RUN_TYPE} sensitivity run ${RUN_ID}."
echo "Results will be written to ${RUN_DIR}."

Rscript run_sensitivity.R \
  outdir="${RUN_DIR}" \
  mode="${RUN_TYPE}" \
  run_id="${RUN_ID}" \
  reps="${REPS}" \
  n_vec="${N_VEC}" \
  N_vec="${OBS_VEC}" \
  m_vec="${M_VEC}" \
  kernels=gauss,epan,quar \
  bw_multipliers=0.5,0.75,1,1.25,1.5 \
  cores="${SLURM_CPUS_PER_TASK}" \
  seed_base=2000 \
  setting=1 \
  decayR=1 \
  latent=FALSE \
  methodXi=CE \
  t_design=random \
  t_grid_len=51 \
  eps_u=0.05 \
  u_cap=0.95 \
  Z_cap=8 \
  beta_a=0.6,3 \
  beta_b=0.6,3 \
  mod_script="${MODULE_SCRIPT}"

if [[ ! -s "${RUN_DIR}/raw_results.rds" || ! -s "${RUN_DIR}/RUN_COMPLETE" ]]; then
  echo "ERROR: the simulation did not create complete combined output." >&2
  exit 5
fi
if [[ ! -d "${RUN_DIR}/replications" ]]; then
  echo "ERROR: the simulation did not create the replication directory." >&2
  exit 5
fi
shopt -s nullglob
replication_files=("${RUN_DIR}"/replications/rep_[0-9][0-9][0-9].rds)
shopt -u nullglob
if (( ${#replication_files[@]} != REPS )); then
  echo "ERROR: expected ${REPS} replication files; found ${#replication_files[@]}." >&2
  exit 5
fi

Rscript aggregate_results.R \
  input_dir="${RUN_DIR}" \
  outdir="${RUN_DIR}/aggregated" \
  mode="${RUN_TYPE}" \
  strict=TRUE \
  write_table="${WRITE_TABLE}" \
  table_out="${RUN_DIR}/aggregated/sensitivity_table.tex"

if grep -qx 'status=COMPLETE_WITH_FIT_FAILURES' "${RUN_DIR}/aggregated/numerical_status.txt"; then
  COMPLETION_KIND="FIT_FAILURES"
elif ! grep -qx 'status=COMPLETE' "${RUN_DIR}/aggregated/numerical_status.txt"; then
  echo "ERROR: aggregation produced an unrecognized numerical status." >&2
  exit 6
fi

echo "Sensitivity run ${RUN_ID} passed strict structural validation."
