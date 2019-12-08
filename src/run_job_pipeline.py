import subprocess


def run_jobs(parts=100,
             amount=637):
    # 637 drugs in intersect
    # 186135 for all drugs in STITCH -> set human_only in msa to wrong
    # 1430 drugs in SIDER
    slurm_path = "../SLURM_JOBS/"

    sep_intervals = [(int(amount / parts * i), int(amount / parts * (i + 1))) for i in range(parts)]

    for start, end in sep_intervals:
        preface_script = '''#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J MSA
#SBATCH -o MultiSequenceAlignment.%J.out
#SBATCH -e MultiSequenceAlignment.%J.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=20

#run the application:
module load anaconda3/4.4.0
module load gcc/6.4.0
source /home/${USER}/.bashrc
conda activate ~/.conda/envs/dti/
        
'''
        filename = slurm_path+str(start)+"_jobscript.sh"
        with open(file=filename, mode='w') as f:
            f.write(preface_script)
            f.write("python3 alignment_pipeline.py " + str(start) + " " + str(end))

        subprocess.call("sbatch "+filename, shell=True)


def cancel_jobs():
    subprocess.call("squeue | grep 'hinnertr' > cancellist", shell=True)
    filename = "cancellist"
    with open(file=filename, mode='r') as f:
        for line in f:
            split_line = line.split(' ')
            for ele in split_line:
                if ele!='':
                    print(ele)
                    subprocess.call("scancel " + ele, shell=True)
                    break

def submit_jobscript_n_times(n):
    for _ in range(n):
        subprocess.call("sbatch jobscript.sh", shell=True)

if __name__ == '__main__':
    # run_jobs(parts=60, amount=646)
    cancel_jobs()
    submit_jobscript_n_times(4)

    pass