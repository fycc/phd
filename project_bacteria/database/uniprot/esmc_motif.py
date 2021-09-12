#!/lustre03/project/6015106/fyc/venv/bin/python
#SBATCH --account=def-yxia
#SBATCH --time=6:00:00
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=8000M
#SBATCH --job-name=esmc_motif
#SBATCH --output=%x-%j.out

from Bio import SeqIO
import os, subprocess
import re

path = '/lustre03/project/6015106/fyc/database/uniprot'

os.chdir(path)
os.environ['IPATH_NO_CPUAFFINITY'] = '1'
os.environ['OMP_NUM_THREADS'] = '16'

motifs = [re.compile(m) for m in open('/lustre03/project/6015106/fyc/database/3did_elm/all_motifs.txt').read().splitlines()]

def elm(infile):
	record_iterator = SeqIO.parse(infile, 'fasta')
	with open(infile.replace('.fasta', '_motif.csv'), 'a') as outfile:
		for record in record_iterator:
			mm = []
			for motif in motifs:
				mm.append(str(len(re.findall(motif, str(record.seq)))))
			outfile.write(','.join([record.id.split('|')[2], str(len(record.seq)), ";".join(mm)]) + "\n")

elm('esmc.fasta')
