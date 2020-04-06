import subprocess


def run_jobs(parts=100,
             amount=637):
    # 637 drugs in intersect
    # 186135 for all drugs in STITCH -> set human_only in msa to wrong
    # 1430 drugs in SIDER
    # 2254 union SIDER and merged graph
    slurm_path = "../SLURM_JOBS/"

    sep_intervals = [(int(amount / parts * i), int(amount / parts * (i + 1))) for i in range(parts)]

    for start, end in sep_intervals:
        preface_script = '''#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J MAFFTMSA
#SBATCH -o MAFFTMultiSequenceAlignment.%J.out
#SBATCH -e MAFFTMultiSequenceAlignment.%J.err
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


def submit_gpu_job(num_proteins=-1,
                   epochs=None,
                   batch_size=None,
                   arch='SimpleGCN',
                   node_features='simple'):
    jobname = arch+'_'+node_features
    preface_script = '''#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J {jobname}
#SBATCH -o jobscript_outputs/{jobname}.%J.out
#SBATCH -e jobscript_outputs/{jobname}.%J.err
#SBATCH --time=3-00:00:00
#SBATCH --gres=gpu:v100:4
#SBATCH --mem=300G
#SBATCH --constraint=[gpu]

#run the application:
module load anaconda3/4.4.0
source /home/hinnertr/.bashrc
conda activate ~/.conda/envs/dti/

module load cuda/10.0.130

export CUDA_VISIBLE_DEVICES=0,1,2,3
python3 torch_dti_predictor.py '''.format(jobname=jobname)
    preface_script += "--num_proteins {num_prots} ".format(num_prots=str(num_proteins))
    if epochs:
        preface_script += "--num_epochs={num_epochs} ".format(num_epochs=str(epochs))
    if batch_size:
        preface_script += "--batch_size={batch_size} ".format(batch_size=str(batch_size))
    preface_script += "--num_folds 5 "

    filename = '../SLURM_jobs/'+jobname+'_jobscript.sh'
    with open(file=filename, mode='w') as f:
        f.write(preface_script)
    subprocess.call("sbatch " + filename, shell=True)


if __name__ == '__main__':
    # run_jobs(parts=150, amount=2254)
    # cancel_jobs()
    # submit_jobscript_n_times(50)

    submit_gpu_job(epochs=30, batch_size=32)
    submit_gpu_job(num_proteins=4000, epochs=30, batch_size=128)
