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
                   arch='GCNConv',
                   days=2,
                   num_gpus=4,
                   node_features='MolPred',
                   fold=-1,
                   mode='standard'):
    jobname = arch+'_'+node_features+'_'+mode
    preface_script = '''#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J {jobname}
#SBATCH -o jobscript_outputs/{jobname}.%J.out
#SBATCH -e jobscript_outputs/{jobname}.%J.err
#SBATCH --time={days}-00:00:00
#SBATCH --gres=gpu:v100:{num_gpus}
#SBATCH --mem=500G
#SBATCH --constraint=[gpu]
# SBATCH --sockets-per-node=1
# SBATCH --gpus-per-socket={num_gpus}
#SBATCH --cpus-per-gpu=6

#run the application:
module load anaconda3/4.4.0
source /home/hinnertr/.bashrc
conda activate ~/.conda/envs/dti/

module load cuda/10.0.130

export CUDA_VISIBLE_DEVICES={vis_dev}
python3 torch_dti_predictor.py '''.format(jobname=jobname, days=str(days), num_gpus=str(num_gpus), vis_dev=str(list(range(num_gpus)))[1:-1].replace(' ', ''))
    preface_script += "--num_proteins {num_prots} ".format(num_prots=str(num_proteins))
    if epochs:
        preface_script += "--num_epochs={num_epochs} ".format(num_epochs=str(epochs))
    if batch_size:
        preface_script += "--batch_size={batch_size} ".format(batch_size=str(batch_size))

    preface_script += "--fold {fold} ".format(fold=str(fold))
    preface_script += "--arch {arch} ".format(arch=arch)
    preface_script += "--mode {mode} ".format(mode=mode)

    filename = '../SLURM_JOBS/'+jobname+'_jobscript.sh'
    with open(file=filename, mode='w') as f:
        f.write(preface_script)
    subprocess.call("sbatch " + filename, shell=True)


def submit_mol_pred_gpu_job(num_proteins=-1,
                            epochs=None,
                            batch_size=None,
                            model='protein',
                            days=1,
                            num_gpus=4,
                            fold=-1):
    jobname = 'Mol_Pred'
    preface_script = '''#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J {jobname}
#SBATCH -o jobscript_outputs/{jobname}.%J.out
#SBATCH -e jobscript_outputs/{jobname}.%J.err
#SBATCH --time={days}-00:00:00
#SBATCH --gres=gpu:v100:{num_gpus}
#SBATCH --mem-per-gpu=60G
#SBATCH --constraint=[gpu]
# SBATCH --sockets-per-node=1
# SBATCH --gpus-per-socket={num_gpus}
#SBATCH --cpus-per-gpu=6

#run the application:
module load anaconda3/4.4.0
source /home/hinnertr/.bashrc
conda activate ~/.conda/envs/dti/

module load cuda/10.0.130

python3 molecular_predictor.py '''.format(jobname=jobname, days=str(days), num_gpus=str(num_gpus))
    preface_script += "--num_proteins {num_prots} ".format(num_prots=str(num_proteins))
    if epochs:
        preface_script += "--num_epochs={num_epochs} ".format(num_epochs=str(epochs))
    if batch_size:
        preface_script += "--batch_size={batch_size} ".format(batch_size=str(batch_size))
    preface_script += "--fold {fold} ".format(fold=str(fold))
    preface_script += "--model {model} ".format(model=model)

    filename = '../SLURM_JOBS/'+jobname+'_jobscript.sh'
    with open(file=filename, mode='w') as f:
        f.write(preface_script)
    subprocess.call("sbatch " + filename, shell=True)


if __name__ == '__main__':
    # run_jobs(parts=150, amount=2254)
    # cancel_jobs()
    # submit_jobscript_n_times(50)

    # cancel_jobs()
    # 'ChebConv','GraphConv', 'TAGConv', 'ARMAConv', 'SGConv', 'FeaStConv'
    for arch in ['GCNConv','SAGEConv', 'GATConv']:
        # submit_gpu_job(epochs=6, batch_size=80, days=2, arch=arch, mode='protein_drughub')
        # submit_gpu_job(epochs=6, batch_size=80, days=2, arch=arch, mode='drug_drughub')

        for fold in range(1, 3):
            # submit_gpu_job(epochs=30, batch_size=32, days=2, arch=arch, mode='drug', fold=fold)

            submit_gpu_job(epochs=10, batch_size=8, days=2, arch=arch, fold=fold, num_gpus=4)
            # submit_gpu_job(num_proteins=4000, epochs=30, batch_size=64, arch=arch)
            # submit_gpu_job(num_proteins=1000, epochs=30, batch_size=256, arch=arch)

            submit_gpu_job(epochs=10, batch_size=8, days=2, arch='Res'+arch, fold=fold, num_gpus=4)
            # submit_gpu_job(num_proteins=4000, epochs=30, batch_size=64, arch='Res'+arch)
            # submit_gpu_job(num_proteins=1000, epochs=30, batch_size=256, arch='Res'+arch)

