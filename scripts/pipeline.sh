#!/bin/bash

# Gaussian Calculation Pipeline Manager
# =======================================
# This script checks for completed Gaussian calculations and submits the next job in the sequence
# Works for all 4-character PDB directories in the current folder
#
# USAGE:
#   ./gaussian_pipeline_manager.sh
#
# CALCULATION ORDER:
#   01 → 01L2 → [05re, 12re] → 05 → 05L4
#                            → 12 → 12L6 → 16re → 16 → 16L8
#
# DIRECTORY STRUCTURE:
#   {pdb}/           (e.g., 1a2f/)
#     ├── {pdb}01/   Contains: 01, 01L2 calculations
#     ├── {pdb}05/   Contains: 05re, 05, 05L4 calculations
#     ├── {pdb}12/   Contains: 12re, 12, 12L6 calculations
#     └── {pdb}16/   Contains: 16re, 16, 16L8 calculations
#
# COMPLETION CRITERIA:
#   1. Checkpoint file exists and is larger than 20MB
#   2. Log file contains "Normal termination of Gaussian" in last 10 lines
#
# NOTES:
#   - States 05 and 12 run in parallel after 01L2 completes
#   - State 16 begins only after 12L6 completes
#   - Checkpoint files are automatically copied between state folders as needed
#
# CHECKPOINT COPIES (automatic):
#   - 01L2 → 05re: Copies final .chk from {pdb}01/ to {pdb}05/
#   - 01L2 → 12re: Copies final .chk from {pdb}01/ to {pdb}12/
#   - 12L6 → 16re: Copies final .chk from {pdb}12/ to {pdb}16/

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Minimum checkpoint file size (20MB in bytes)
MIN_CHK_SIZE=$((20 * 1024 * 1024))

#################################################
# Function: check_calculation_complete
# Checks if a Gaussian calculation completed successfully
# Args: $1 = subdirectory path, $2 = calculation name (e.g., "01", "05re")
# Returns: 0 if complete, 1 if not complete
#################################################
check_calculation_complete() {
    local subdir="$1"
    local calc_name="$2"
    local pdb=$(basename $(dirname "$subdir"))

    # Construct file paths
    local com_file="${subdir}/${pdb}${calc_name}.com"
    local log_file="${subdir}/${pdb}${calc_name}.log"

    # Check if .com file exists
    if [[ ! -f "$com_file" ]]; then
        echo -e "${YELLOW}    [SKIP] ${com_file} not found${NC}"
        return 1
    fi

    # Extract checkpoint filename from .com file
    local chk_file=$(grep "^%chk=" "$com_file" | tail -1 | cut -d'=' -f2)
    if [[ -z "$chk_file" ]]; then
        echo -e "${RED}    [ERROR] No checkpoint file specified in ${com_file}${NC}"
        return 1
    fi

    local chk_path="${subdir}/${chk_file}"

    # Check 1: Does checkpoint file exist and is it > 20MB?
    if [[ ! -f "$chk_path" ]]; then
        echo -e "${YELLOW}    [WAITING] Checkpoint ${chk_file} not found${NC}"
        return 1
    fi

    local chk_size=$(stat -f%z "$chk_path" 2>/dev/null || stat -c%s "$chk_path" 2>/dev/null)
    if [[ $chk_size -lt $MIN_CHK_SIZE ]]; then
        echo -e "${YELLOW}    [WAITING] Checkpoint ${chk_file} too small ($(($chk_size / 1024 / 1024))MB < 20MB)${NC}"
        return 1
    fi

    # Check 2: Does log file contain "Normal termination" in last 10 lines?
    if [[ ! -f "$log_file" ]]; then
        echo -e "${YELLOW}    [WAITING] Log file ${pdb}${calc_name}.log not found${NC}"
        return 1
    fi

    if tail -10 "$log_file" | grep -q "Normal termination of Gaussian"; then
        echo -e "${GREEN}    [COMPLETE] ${calc_name} finished successfully${NC}"
        return 0
    else
        echo -e "${RED}    [FAILED] ${calc_name} did not terminate normally${NC}"
        return 1
    fi
}

#################################################
# Function: copy_checkpoint_file
# Copies a checkpoint file from source to destination directory
# Args: $1 = source subdirectory, $2 = source calc name,
#       $3 = destination subdirectory
# Returns: 0 if successful, 1 if failed
#################################################
copy_checkpoint_file() {
    local src_subdir="$1"
    local src_calc="$2"
    local dst_subdir="$3"
    local pdb=$(basename $(dirname "$src_subdir"))

    # Get the checkpoint file name from source calculation
    local src_com="${src_subdir}/${pdb}${src_calc}.com"
    local chk_file=$(grep "^%chk=" "$src_com" | tail -1 | cut -d'=' -f2)

    if [[ -z "$chk_file" ]]; then
        echo -e "${RED}    [ERROR] Cannot find checkpoint file in ${src_com}${NC}"
        return 1
    fi

    local src_chk="${src_subdir}/${chk_file}"
    local dst_chk="${dst_subdir}/${chk_file}"

    # Check if source checkpoint exists
    if [[ ! -f "$src_chk" ]]; then
        echo -e "${RED}    [ERROR] Source checkpoint ${src_chk} not found${NC}"
        return 1
    fi

    # Copy if destination doesn't exist or is different
    if [[ ! -f "$dst_chk" ]] || ! cmp -s "$src_chk" "$dst_chk"; then
        echo -e "${BLUE}    [COPY] ${chk_file} from ${src_subdir} to ${dst_subdir}${NC}"
        cp "$src_chk" "$dst_chk"
        if [[ $? -eq 0 ]]; then
            echo -e "${GREEN}    [COPY OK] Checkpoint copied successfully${NC}"
            return 0
        else
            echo -e "${RED}    [COPY FAILED] Failed to copy checkpoint${NC}"
            return 1
        fi
    else
        echo -e "${GREEN}    [COPY OK] Checkpoint already present in destination${NC}"
        return 0
    fi
}

#################################################
# Function: is_job_running
# Checks if a job is currently running in SLURM queue
# Args: $1 = subdirectory path, $2 = calculation name
# Returns: 0 if running, 1 if not running
#################################################
is_job_running() {
    local subdir="$1"
    local calc_name="$2"
    local pdb=$(basename $(dirname "$subdir"))

    # Look for the job in the SLURM queue
    # The job script name would be g16_{pdb}{calc_name}.sh
    local job_pattern="g16_${pdb}${calc_name}"

    # Check if job is in queue (running or pending)
    if squeue -u $USER -o "%.100j" 2>/dev/null | grep -q "${job_pattern}"; then
        return 0  # Job is running
    else
        return 1  # Job is not running
    fi
}

#################################################
# Function: submit_calculation
# Submits a Gaussian calculation using sbatch
# Args: $1 = subdirectory path, $2 = calculation name (e.g., "01", "05re")
#################################################
submit_calculation() {
    local subdir="$1"
    local calc_name="$2"
    local pdb=$(basename $(dirname "$subdir"))

    local submit_script="${subdir}/g16_${pdb}${calc_name}.sh"

    if [[ ! -f "$submit_script" ]]; then
        echo -e "${RED}    [ERROR] Submission script ${submit_script} not found${NC}"
        return 1
    fi

    # Check if job is already running in SLURM queue
    if is_job_running "$subdir" "$calc_name"; then
        echo -e "${YELLOW}    [RUNNING] ${calc_name} is already in SLURM queue${NC}"
        return 1
    fi

    echo -e "${BLUE}    [SUBMIT] Submitting ${calc_name} via sbatch${NC}"
    cd "$subdir" || return 1
    sbatch "g16_${pdb}${calc_name}.sh"
    cd - > /dev/null || return 1
    return 0
}

#################################################
# Function: process_subdirectory
# Processes a specific calculation subdirectory
# Args: $1 = subdirectory path, $2 = subdirectory type (01, 05, 12, 16)
#################################################
process_subdirectory() {
    local subdir="$1"
    local subdir_type="$2"
    local pdb=$(basename $(dirname "$subdir"))

    echo -e "${BLUE}  Processing ${subdir_type}:${NC}"

    case "$subdir_type" in
        "01")
            # Sequence: 01 -> 01L2
            if check_calculation_complete "$subdir" "01"; then
                if ! check_calculation_complete "$subdir" "01L2"; then
                    submit_calculation "$subdir" "01L2"
                fi
            else
                # 01 not complete, try to submit it
                submit_calculation "$subdir" "01"
            fi
            ;;

        "05")
            # Sequence: 05re -> 05 -> 05L4
            # 05re depends on 01L2 from subdirectory 01
            local subdir_01="${subdir/05/01}"

            if check_calculation_complete "$subdir" "05re"; then
                if check_calculation_complete "$subdir" "05"; then
                    if ! check_calculation_complete "$subdir" "05L4"; then
                        submit_calculation "$subdir" "05L4"
                    fi
                else
                    submit_calculation "$subdir" "05"
                fi
            else
                # Check if 01L2 is complete before starting 05re
                if check_calculation_complete "$subdir_01" "01L2" > /dev/null 2>&1; then
                    # Copy checkpoint from 01L2 to 05 directory
                    if copy_checkpoint_file "$subdir_01" "01L2" "$subdir"; then
                        submit_calculation "$subdir" "05re"
                    else
                        echo -e "${RED}    [ERROR] Failed to copy checkpoint, cannot start 05re${NC}"
                    fi
                else
                    echo -e "${YELLOW}    [WAITING] Need 01L2 to complete before starting 05re${NC}"
                fi
            fi
            ;;

        "12")
            # Sequence: 12re -> 12 -> 12L6
            # 12re depends on 01L2 from subdirectory 01
            local subdir_01="${subdir/12/01}"

            if check_calculation_complete "$subdir" "12re"; then
                if check_calculation_complete "$subdir" "12"; then
                    if ! check_calculation_complete "$subdir" "12L6"; then
                        submit_calculation "$subdir" "12L6"
                    fi
                else
                    submit_calculation "$subdir" "12"
                fi
            else
                # Check if 01L2 is complete before starting 12re
                if check_calculation_complete "$subdir_01" "01L2" > /dev/null 2>&1; then
                    # Copy checkpoint from 01L2 to 12 directory
                    if copy_checkpoint_file "$subdir_01" "01L2" "$subdir"; then
                        submit_calculation "$subdir" "12re"
                    else
                        echo -e "${RED}    [ERROR] Failed to copy checkpoint, cannot start 12re${NC}"
                    fi
                else
                    echo -e "${YELLOW}    [WAITING] Need 01L2 to complete before starting 12re${NC}"
                fi
            fi
            ;;

        "16")
            # Sequence: 16re -> 16 -> 16L8
            # 16re depends on 12L6 from subdirectory 12
            local subdir_12="${subdir/16/12}"

            if check_calculation_complete "$subdir" "16re"; then
                if check_calculation_complete "$subdir" "16"; then
                    if ! check_calculation_complete "$subdir" "16L8"; then
                        submit_calculation "$subdir" "16L8"
                    fi
                else
                    submit_calculation "$subdir" "16"
                fi
            else
                # Check if 12L6 is complete before starting 16re
                if check_calculation_complete "$subdir_12" "12L6" > /dev/null 2>&1; then
                    # Copy checkpoint from 12L6 to 16 directory
                    if copy_checkpoint_file "$subdir_12" "12L6" "$subdir"; then
                        submit_calculation "$subdir" "16re"
                    else
                        echo -e "${RED}    [ERROR] Failed to copy checkpoint, cannot start 16re${NC}"
                    fi
                else
                    echo -e "${YELLOW}    [WAITING] Need 12L6 to complete before starting 16re${NC}"
                fi
            fi
            ;;
    esac
}

#################################################
# Main script
#################################################
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Gaussian Pipeline Manager${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# Find all 4-character directories in current folder
for pdb_dir in */; do
    pdb=$(basename "$pdb_dir")

    # Check if directory name is exactly 4 characters
    if [[ ${#pdb} -ne 4 ]]; then
        continue
    fi

    echo -e "${GREEN}Processing PDB: ${pdb}${NC}"

    # Process each subdirectory in the sequence
    for subdir_type in "01" "05" "12" "16"; do
        subdir="${pdb_dir}${pdb}${subdir_type}"

        if [[ -d "$subdir" ]]; then
            process_subdirectory "$subdir" "$subdir_type"
        else
            echo -e "${YELLOW}  Subdirectory ${subdir} not found, skipping${NC}"
        fi
    done

    echo ""
done

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Pipeline check complete!${NC}"
echo -e "${GREEN}========================================${NC}"

