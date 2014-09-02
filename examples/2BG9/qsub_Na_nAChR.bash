#$ -N BPNa_nAChR
#$ -S /bin/bash
#$ -l mem_free=1024M
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1-99

# :Queuing System: GE
# :Host: ASU Becksteinlab Workstations

# *********** WARNING **********
# THIS SCRIPT WILL ONLY RUN ON THE WORKSTATION QUEUE IN THE BECKSTEINLAB
# You can use it as an example for how to set up the calculation
# on a queuing system with a heterogenous computing environment.
# *********** WARNING **********


# This file is interpolated by Python (i.e. anything with 
# a single or double percentage sign ('%') will be expanded. See
# http://docs.python.org/library/stdtypes.html#string-formatting-operations
# Tip: Do not add or delete percentage signs unless you know 
#      exactly what you are doing.

# ===================
# Prepare binaries
# ===================
# Set platform and architecture

# Init modules
. /nfs/packages/modules/init

# Set platform and architecture
MACHINE_TYPE="`uname -s`_`uname -m`"

# Print some useful info
echo "------------------------------------------------------------"
echo "Execution host:     ${HOSTNAME}"
echo "Date:               $(date)"
echo "User:               ${USER}"
echo "Machine type:       ${MACHINE_TYPE}"
echo "SGE temp directory: ${TMPDIR}" 
echo "Working directory:  ${SGE_O_WORKDIR}" 
echo "------------------------------------------------------------"

# Load modules
case "$MACHINE_TYPE" in
'Linux_x86_64')
	echo "Loading 64-bit Linux modules"
	module load apbs/1.4/64
	module load APBSmem/1.12/64
;; 
'Darwin_i386')
	echo "Loading 32-bit Darwin modules"
	#module load apbs/1.3svn1625/univ
	#module load APBSmem/1.05/univ
	echo "ERROR: No binaries for platform and architecture."
	exit 1
;;
'Darwin_x86_64')
	echo "Loading 64-bit Darwin modules"
	#module load apbs/1.3svn1625/univ
	#module load APBSmem/1.05/univ
	echo "ERROR: No binaries for platform and architecture."
	exit 1
;;
*)
	echo "ERROR: Unknown platform and architecture."
	exit 1
;;
esac

# random delay of 0..MAXDELAY s to avoid nfs storms
: ${MAXDELAY:=60}
delay=$(($RANDOM*${MAXDELAY}/32767))
echo    "Array job number <job>${JOB_ID}.${SGE_TASK_ID}</job> on <host>${HOSTNAME}</host>:"
echo -n "Waiting for ${delay} seconds before accessing a nfs disk..."
sleep ${delay}
echo "done"

# environment variables are used in the downstream scripts to set executables
export APBS=`which apbs`
export DRAW_MEMBRANE2A=`which draw_membrane2a`

# sanity checks
check_passed=True
if [ -n "${APBS}" ] && [ -x "${APBS}" ]; then
   echo "Check passed: Using apbs from ${APBS}"
else
   echo "ERROR: could not find 'apbs' executable."
   check_passed=False
fi
if [ -n "${DRAW_MEMBRANE2A}" ] && [ -x "${DRAW_MEMBRANE2A}" ]; then
   echo "Check passed: Using draw_membrane2a from ${DRAW_MEMBRANE2A}"
else
   echo "ERROR: could not find 'draw_membrane2a' executable."
   check_passed=False
fi
case "${check_passed}" in
     False|false)   exit 2;;
esac

#------------------------------------------------------------
# main
#------------------------------------------------------------

declare -a job

#------------------------------------------------------------
job[1]="Na_nAChR/w0000/job_0000.bash"
job[2]="Na_nAChR/w0001/job_0001.bash"
job[3]="Na_nAChR/w0002/job_0002.bash"
job[4]="Na_nAChR/w0003/job_0003.bash"
job[5]="Na_nAChR/w0004/job_0004.bash"
job[6]="Na_nAChR/w0005/job_0005.bash"
job[7]="Na_nAChR/w0006/job_0006.bash"
job[8]="Na_nAChR/w0007/job_0007.bash"
job[9]="Na_nAChR/w0008/job_0008.bash"
job[10]="Na_nAChR/w0009/job_0009.bash"
job[11]="Na_nAChR/w0010/job_0010.bash"
job[12]="Na_nAChR/w0011/job_0011.bash"
job[13]="Na_nAChR/w0012/job_0012.bash"
job[14]="Na_nAChR/w0013/job_0013.bash"
job[15]="Na_nAChR/w0014/job_0014.bash"
job[16]="Na_nAChR/w0015/job_0015.bash"
job[17]="Na_nAChR/w0016/job_0016.bash"
job[18]="Na_nAChR/w0017/job_0017.bash"
job[19]="Na_nAChR/w0018/job_0018.bash"
job[20]="Na_nAChR/w0019/job_0019.bash"
job[21]="Na_nAChR/w0020/job_0020.bash"
job[22]="Na_nAChR/w0021/job_0021.bash"
job[23]="Na_nAChR/w0022/job_0022.bash"
job[24]="Na_nAChR/w0023/job_0023.bash"
job[25]="Na_nAChR/w0024/job_0024.bash"
job[26]="Na_nAChR/w0025/job_0025.bash"
job[27]="Na_nAChR/w0026/job_0026.bash"
job[28]="Na_nAChR/w0027/job_0027.bash"
job[29]="Na_nAChR/w0028/job_0028.bash"
job[30]="Na_nAChR/w0029/job_0029.bash"
job[31]="Na_nAChR/w0030/job_0030.bash"
job[32]="Na_nAChR/w0031/job_0031.bash"
job[33]="Na_nAChR/w0032/job_0032.bash"
job[34]="Na_nAChR/w0033/job_0033.bash"
job[35]="Na_nAChR/w0034/job_0034.bash"
job[36]="Na_nAChR/w0035/job_0035.bash"
job[37]="Na_nAChR/w0036/job_0036.bash"
job[38]="Na_nAChR/w0037/job_0037.bash"
job[39]="Na_nAChR/w0038/job_0038.bash"
job[40]="Na_nAChR/w0039/job_0039.bash"
job[41]="Na_nAChR/w0040/job_0040.bash"
job[42]="Na_nAChR/w0041/job_0041.bash"
job[43]="Na_nAChR/w0042/job_0042.bash"
job[44]="Na_nAChR/w0043/job_0043.bash"
job[45]="Na_nAChR/w0044/job_0044.bash"
job[46]="Na_nAChR/w0045/job_0045.bash"
job[47]="Na_nAChR/w0046/job_0046.bash"
job[48]="Na_nAChR/w0047/job_0047.bash"
job[49]="Na_nAChR/w0048/job_0048.bash"
job[50]="Na_nAChR/w0049/job_0049.bash"
job[51]="Na_nAChR/w0050/job_0050.bash"
job[52]="Na_nAChR/w0051/job_0051.bash"
job[53]="Na_nAChR/w0052/job_0052.bash"
job[54]="Na_nAChR/w0053/job_0053.bash"
job[55]="Na_nAChR/w0054/job_0054.bash"
job[56]="Na_nAChR/w0055/job_0055.bash"
job[57]="Na_nAChR/w0056/job_0056.bash"
job[58]="Na_nAChR/w0057/job_0057.bash"
job[59]="Na_nAChR/w0058/job_0058.bash"
job[60]="Na_nAChR/w0059/job_0059.bash"
job[61]="Na_nAChR/w0060/job_0060.bash"
job[62]="Na_nAChR/w0061/job_0061.bash"
job[63]="Na_nAChR/w0062/job_0062.bash"
job[64]="Na_nAChR/w0063/job_0063.bash"
job[65]="Na_nAChR/w0064/job_0064.bash"
job[66]="Na_nAChR/w0065/job_0065.bash"
job[67]="Na_nAChR/w0066/job_0066.bash"
job[68]="Na_nAChR/w0067/job_0067.bash"
job[69]="Na_nAChR/w0068/job_0068.bash"
job[70]="Na_nAChR/w0069/job_0069.bash"
job[71]="Na_nAChR/w0070/job_0070.bash"
job[72]="Na_nAChR/w0071/job_0071.bash"
job[73]="Na_nAChR/w0072/job_0072.bash"
job[74]="Na_nAChR/w0073/job_0073.bash"
job[75]="Na_nAChR/w0074/job_0074.bash"
job[76]="Na_nAChR/w0075/job_0075.bash"
job[77]="Na_nAChR/w0076/job_0076.bash"
job[78]="Na_nAChR/w0077/job_0077.bash"
job[79]="Na_nAChR/w0078/job_0078.bash"
job[80]="Na_nAChR/w0079/job_0079.bash"
job[81]="Na_nAChR/w0080/job_0080.bash"
job[82]="Na_nAChR/w0081/job_0081.bash"
job[83]="Na_nAChR/w0082/job_0082.bash"
job[84]="Na_nAChR/w0083/job_0083.bash"
job[85]="Na_nAChR/w0084/job_0084.bash"
job[86]="Na_nAChR/w0085/job_0085.bash"
job[87]="Na_nAChR/w0086/job_0086.bash"
job[88]="Na_nAChR/w0087/job_0087.bash"
job[89]="Na_nAChR/w0088/job_0088.bash"
job[90]="Na_nAChR/w0089/job_0089.bash"
job[91]="Na_nAChR/w0090/job_0090.bash"
job[92]="Na_nAChR/w0091/job_0091.bash"
job[93]="Na_nAChR/w0092/job_0092.bash"
job[94]="Na_nAChR/w0093/job_0093.bash"
job[95]="Na_nAChR/w0094/job_0094.bash"
job[96]="Na_nAChR/w0095/job_0095.bash"
job[97]="Na_nAChR/w0096/job_0096.bash"
job[98]="Na_nAChR/w0097/job_0097.bash"
job[99]="Na_nAChR/w0098/job_0098.bash" 
#------------------------------------------------------------

run_d=$(dirname ${job[${SGE_TASK_ID}]})
script=$(basename ${job[${SGE_TASK_ID}]})
 
cd ${run_d} || { echo "Failed to cd ${run_d}. Abort."; exit 1; }
. ./${script}
