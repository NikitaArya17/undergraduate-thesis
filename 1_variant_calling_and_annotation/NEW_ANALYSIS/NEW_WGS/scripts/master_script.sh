#!/bin/bash

# To submit 5 jobs sequentially
SCRIPT_DIR="/home/nikita.arya/GG20_WGS_NEW/scripts"
SCRIPT1="${SCRIPT_DIR}/1_fastqc_reads.sh"
SCRIPT2="${SCRIPT_DIR}/2_bwa_alignment.sh"
SCRIPT3="${SCRIPT_DIR}/3_mark_duplicates.sh"
SCRIPT4="${SCRIPT_DIR}/4_samtools_stats.sh"
SCRIPT5="${SCRIPT_DIR}/5_variant_calling.sh"
JOB_LOG="${SCRIPT_DIR}/job_log.txt"

JID1=$(sbatch --parsable ${SCRIPT1})
echo "Submitted Step 1 with Job ID: $JID1" > ${JOB_LOG}

JID2=$(sbatch --parsable --dependency=afterok:$JID1 ${SCRIPT2})
echo "Submitted Step 2 with Job ID: $JID2 (waiting for $JID1)" >> ${JOB_LOG}

JID3=$(sbatch --parsable --dependency=afterok:$JID2 ${SCRIPT3})
echo "Submitted Step 3 with Job ID: $JID3 (waiting for $JID2)" >> ${JOB_LOG}

JID4=$(sbatch --parsable --dependency=afterok:$JID3 ${SCRIPT4})
echo "Submitted Step 4 with Job ID: $JID4 (waiting for $JID3)" >> ${JOB_LOG}

JID5=$(sbatch --parsable --dependency=afterok:$JID4 ${SCRIPT5})
echo "Submitted Step 5 with Job ID: $JID5 (waiting for $JID4)" >> ${JOB_LOG}
