import os
import subprocess

program = '/home/ramon/MEGAsync/prodgraph/code/RandomWalk/Debug/RandomWalk'
output_folder = '/home/ramon/MEGAsync/prodgraph/code/RandomWalk/Debug/'

filename= '/home/ramon/MEGAsync/prodgraph/code/RandomWalk/Debug/small.csv'

#define RANDOM_ALGORITHM 1
#define GORI_ALGORITHM 2

#define CLIQUE_STRATEGY 1
#define INTERSECTION_STRATEGY 2

k_fold = 5

def run_experiments(k_fold, algorithm, strategy):
    output_file = output_folder+str(algorithm)+'_'+str(strategy)
    f = open(output_file, 'w')
    command = [program, filename, k_fold, algorithm, strategy]
    subprocess.call(command, stdout=f)
    f.close()

run_experiments(k_fold, 1, 1)
run_experiments(k_fold, 1, 2)
run_experiments(k_fold, 2, 1)
