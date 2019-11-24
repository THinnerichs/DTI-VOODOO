


def run_jobs(parts=100,
             amount=186135):

    sep_intervals = [(int(amount / parts * i), int(amount / parts * (i + 1))) for i in range(parts)]

    for start, end in sep_intervals:
        f.write("srun python3 alignment_pipeline.py " + str(start) + " " + str(end) + '\n')
