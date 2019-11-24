import subprocess


def run_jobs(parts=100,
             amount=186135):
    slurm_path = "../SLURM_JOBS/"

    sep_intervals = [(int(amount / parts * i), int(amount / parts * (i + 1))) for i in range(parts)]

    for start, end in sep_intervals:
        preface_script = '''#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J MultiSequenceAlignment
#SBATCH -o MultiSequenceAlignment.%J.out
#SBATCH -e MultiSequenceAlignment.%J.err
#SBATCH --time=72:00:00
#SBATCH --mem=240G
#SBATCH --constraint=[intel]
#SBATCH --cpus-per-task=40

#run the application:
module load anaconda3/4.4.0
source /home/${USER}/.bashrc
conda activate ~/.conda/envs/dti/
        
'''
        filename = slurm_path+str(start)+"_jobscript.sh"
        with open(file=filename, mode='w') as f:
            f.write(preface_script)
            f.write("python3 alignment_pipeline.py " + str(start) + " " + str(end))

        subprocess.call("sbatch "+filename, shell=True)


def cancel_jobs():
    filename = "cancellist"
    with open(file=filename, mode='r') as f:
        for line in f:
            split_line = line.split('\t')
            print(split_line)
            # subprocess.call("scancel "+split_line[0])
if __name__ == '__main__':
    # run_jobs()
    cancel_jobs()