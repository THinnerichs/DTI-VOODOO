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


def cancel_certain_jobs():
    subprocess.call("squeue | grep 'hinnertr' > cancellist", shell=True)
    filename = "cancellist"
    with open(file=filename, mode='r') as f:
        for line in f:
            split_line = line.split(' ')
            job_num = None
            cancel = False
            first = True
            for ele in split_line:
                if ele!='':
                    if first:
                        job_num = ele
                        first = False
                    if ele in 'ChebConv' or ele in 'GraphConv' or ele in 'TAGConv' or ele in 'ARMAConv' or ele in 'SGConv' or ele in 'FeaStConv' or\
                        ele in  'ResChebConv' or ele in 'ResGraphConv' or ele in 'ResTAGConv' or ele in 'ResARMAConv' or ele in 'ResSGConv' or ele in 'ResFeaStConv':
                        cancel = True
            if cancel:
                print(job_num)
                subprocess.call("scancel " + job_num, shell=True)

def cancel_jobs():
    subprocess.call("squeue | grep 'hinnertr' > cancellist", shell=True)
    filename = "cancellist"
    with open(file=filename, mode='r') as f:
        for line in f:
            split_line = line.split(' ')
            for ele in split_line:
                if ele != '':
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
                   mem=300,
                   neg_sample_ratio=1.0,
                   mode='standard',
                   pretrain=True,
                   model_id=''):
    jobname = arch+'_'+node_features+'_'+mode
    preface_script = '''#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J {jobname}
#SBATCH -o jobscript_outputs/{jobname}.%J.out
#SBATCH -e jobscript_outputs/{jobname}.%J.err
#SBATCH --time={days}-00:00:00
#SBATCH --gres=gpu:v100:{num_gpus}
#SBATCH --mem={mem}G
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
python3 torch_dti_predictor.py '''.format(jobname=jobname,
                                          days=str(days),
                                          num_gpus=str(num_gpus),
                                          vis_dev=str(list(range(num_gpus)))[1:-1].replace(' ', ''),
                                          mem=str(mem))
    preface_script += "--num_proteins {num_prots} ".format(num_prots=str(num_proteins))
    if epochs:
        preface_script += "--num_epochs={num_epochs} ".format(num_epochs=str(epochs))
    if batch_size:
        preface_script += "--batch_size={batch_size} ".format(batch_size=str(batch_size))

    preface_script += "--fold {fold} ".format(fold=str(fold))
    preface_script += "--arch {arch} ".format(arch=arch)
    preface_script += "--mode {mode} ".format(mode=mode)
    preface_script += "--neg_sample_ratio {neg_sample_ratio} ".format(neg_sample_ratio=str(neg_sample_ratio))
    preface_script += "--pretrain {pretrain} ".format(pretrain=pretrain)
    preface_script += "--model_id {model_id} ".format(model_id=model_id)

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

def submit_protfunc_pred_job(num_proteins=-1,
                            epochs=None,
                            batch_size=None,
                            model='protein',
                            days=1,
                            num_gpus=4,
                            fold=-1):
    jobname = 'ProFunc'
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

python3 protein_function_predictor.py '''.format(jobname=jobname, days=str(days), num_gpus=str(num_gpus))
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

    # for fold in range(1,6):
        # submit_protfunc_pred_job(epochs=30, batch_size=131072, fold=fold, num_gpus=4, days=1)

    # 'ChebConv','GraphConv', 'TAGConv', 'ARMAConv', 'SGConv', 'FeaStConv', 'SAGEConv', 'GATConv
    for arch in ['GCNConv']:
    # for arch in ['ChebConv','GraphConv', 'TAGConv', 'ARMAConv', 'SGConv', 'FeaStConv']:
        # submit_gpu_job(epochs=20, batch_size=160, mem=360, days=1, arch=arch, mode='protein_drughub', num_gpus=4, neg_sample_ratio=0.05)
        # submit_gpu_job(epochs=20, batch_size=160, mem=360, days=1, arch='Res'+arch, mode='protein_drughub', num_gpus=4, neg_sample_ratio=0.05)
        # submit_gpu_job(epochs=20, batch_size=160, mem=360, days=1, arch=arch, mode='drug_drughub', num_gpus=4, neg_sample_ratio=0.05)
        # submit_gpu_job(epochs=20, batch_size=160, mem=360, days=1, arch='Res'+arch, mode='drug_drughub', num_gpus=4, neg_sample_ratio=0.05)

        for fold in range(1, 3):
            # submit_gpu_job(epochs=30, batch_size=32, days=2, arch=arch, mode='drug', fold=fold, num_gpus=2, neg_sample_ratio=0.05)
            # submit_gpu_job(epochs=30, batch_size=32, days=2, arch='Res'+arch, mode='drug', fold=fold, num_gpus=2, neg_sample_ratio=0.05)

            submit_gpu_job(epochs=18, batch_size=160, mem=360, days=1, arch=arch, fold=fold, num_gpus=4, neg_sample_ratio=0.05, model_id='Shallow', node_features='ProtFunc')
            # submit_gpu_job(epochs=30, batch_size=160, mem=360, days=2, arch=arch, fold=fold, num_gpus=4, neg_sample_ratio=0.1)
            # submit_gpu_job(num_proteins=4000, epochs=30, batch_size=64, arch=arch)
            # submit_gpu_job(num_proteins=1000, epochs=30, batch_size=256, arch=arch)

            # submit_gpu_job(epochs=10, batch_size=140, mem=360, days=1, arch='Res'+arch, fold=fold, num_gpus=4, neg_sample_ratio=0.05, model_id='ResBlocks')
            # submit_gpu_job(epochs=30, batch_size=160, mem=360, days=2, arch='Res'+arch, fold=fold, num_gpus=4, neg_sample_ratio=0.1)
            # submit_gpu_job(num_proteins=4000, epochs=30, batch_size=64, arch='Res'+arch)
            # submit_gpu_job(num_proteins=1000, epochs=30, batch_size=256, arch='Res'+arch)
