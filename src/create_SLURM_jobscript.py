import numpy as np


def distribute(parts=100,
               amount=186135):
    array = np.array(list(range(amount)))

    # sep = np.split(array, parts)

    # sep_intervals = [(a.min(), a.max()) for a in sep]
    sep_intervals = [(int(amount/parts*i),int(amount/parts * (i+1))) for i in range(parts)]

    with open(file="jobscript.sh", mode='a') as f:
        for start, end in sep_intervals:
            f.write("srun python3 alignment_pipeline.py "+str(start) +" "+str(end)+'\n')


if __name__ == '__main__':
    distribute(parts=100)