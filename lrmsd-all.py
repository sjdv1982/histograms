"""
This scripts assume that each docking directory contains:
  [motif].dat : the final dat file
  frag{}r.pdb : the bound form
  [motif].list : the list of fragments fro mthe library (coarse-grain)
"""
import sys, os, subprocess, multiprocessing
docking_dir = sys.argv[1]
if len(sys.argv) > 2:
    nprocesses = int(sys.argv[2])
else:
    nprocesses = -1
os.chdir(docking_dir)
assert os.path.exists("boundfrag.list")
boundfrag = [(l.split()[0], l.split()[1]) for l in open("boundfrag.list") if len(l.strip())]
jobs = []
for frag, motif in boundfrag:
    bound = "frag{}r.pdb".format(frag)
    docked_poses = "{}.dat".format(motif)
    ens_list = "{}.list".format(motif)
    assert os.path.exists(bound), bound
    assert os.path.exists(ens_list), ens_list
    assert os.path.exists(docked_poses), docked_poses
    cmd = "python2 $ATTRACTDIR/lrmsd.py {0} `head -1 {2}` {1} --ens 2 {2} --allatoms > frag{3}.rmsd".format(
        docked_poses, bound, ens_list, frag
    )
    jobs.append(cmd)

def run_job(cmd):
    print(cmd)
    subprocess.run(cmd, check=True, shell=True)

pool = multiprocessing.Pool(nprocesses)
pool.map(run_job, jobs)