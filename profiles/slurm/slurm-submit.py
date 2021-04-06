#!/usr/bin/env python3
import sys
import os
import re
import argparse
import subprocess

from snakemake.utils import read_job_properties

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument(
    "--help", help="Display help message.", action="store_true")
parser.add_argument(
    "positional", action="append",
    nargs="?", metavar="POS",
    help="additional arguments not in slurm parser group to pass to sbatch")

# A subset of SLURM-specific arguments
slurm_parser = parser.add_argument_group("slurm-specific arguments")
slurm_parser.add_argument(
    "-a", "--array", help="job array index values")
slurm_parser.add_argument(
    "-A", "--account", help="charge job to specified account")
slurm_parser.add_argument(
    "--begin", help="defer job until HH:MM MM/DD/YY")
slurm_parser.add_argument(
    "-c", "--cpus-per-task", help="number of cpus required per task")
slurm_parser.add_argument(
    "-d", "--dependency",
    help="defer job until condition on jobid is satisfied")
slurm_parser.add_argument(
    "-D", "--workdir", help="set working directory for batch script")
slurm_parser.add_argument(
    "-J", "--job-name", help="name of job")
slurm_parser.add_argument(
    "--mail-type", help="notify on state change: BEGIN, END, FAIL or ALL")
slurm_parser.add_argument(
    "--mail-user", help="who to send email notification for job state changes")
slurm_parser.add_argument(
    "-n", "--ntasks", help="number of tasks to run")
slurm_parser.add_argument(
    "-N", "--nodes", help="number of nodes on which to run (N = min[-max])")
slurm_parser.add_argument(
    "-p", "--partition", help="partition requested")
slurm_parser.add_argument(
    "-q", "--qos", help="quality of service")
slurm_parser.add_argument(
    "-Q", "--quiet", help="quiet mode (suppress informational messages)")
slurm_parser.add_argument(
    "-t", "--time", help="time limit")
slurm_parser.add_argument(
    "--wrap", help="wrap command string in a sh script and submit")
slurm_parser.add_argument(
    "-C", "--constraint", help="specify a list of constraints")
slurm_parser.add_argument(
    "--mem", help="minimum amount of real memory")

args = parser.parse_args()

if args.help:
    parser.print_help()
    sys.exit(0)

jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

extras = ""
if args.positional:
    for m in args.positional:
        if m is not None:
            extras = extras + " " + m

arg_dict = dict(args.__dict__)


# Process resources
if "resources" in job_properties:
    resources = job_properties["resources"]
    if arg_dict["time"] is None:
        if "time_min" in resources:
            arg_dict["time"] = resources["time_min"]
        elif "walltime" in resources:
            arg_dict["time"] = resources["walltime"]
        elif "runtime" in resources:
            arg_dict["time"] = resources["runtime"]
        else:
            arg_dict["time"] = 30
    if arg_dict["mem"] is None:
        if "mem_mb" in resources:
            arg_dict["mem"] = resources["mem_mb"]
        elif "mem" in resources:
            arg_dict["mem"] = resources["mem"]
        else:
            arg_dict["mem"] = 1024
    if arg_dict["partition"] is None:
        if arg_dict["time"] < 360:
            arg_dict["partition"] = "shortq"
        elif 360 <= arg_dict["time"] < 1440:
            arg_dict["partition"] = "mediumq"
        elif 1440 <= arg_dict["time"] < 10080:
            arg_dict["partition"] = "longq"
        elif 10080 <= arg_dict["time"] < 86400:
            arg_dict["partition"] = "verylongq"
        else:
            raise ValueError(
                "Too much time requested: {}".format(str(arg_dict["time"]))
            )


# Threads
if "threads" in job_properties:
    arg_dict["cpus_per_task"] = job_properties["threads"]

opt_keys = ["array", "account", "begin", "cpus_per_task",
            "dependency", "workdir", "error", "job_name", "mail_type",
            "mail_user", "ntasks", "nodes", "output", "partition",
            "quiet", "time", "wrap", "constraint", "mem"]

arg_dict["output"] = "logs/slurm/slurm-%x-%j-%N.out"
if arg_dict["output"] is not None:
    os.makedirs(os.path.dirname(arg_dict["output"]), exist_ok=True)
arg_dict["error"] = "logs/slurm/slurm-%x-%j-%N.err"
if arg_dict["error"] is not None:
    os.makedirs(os.path.dirname(arg_dict["error"]), exist_ok=True)

arg_dict["mail_type"] = "END,FAIL"
arg_dict["mail_user"] = "{{cookiecutter.mail_user}}"

opts = ""
for k, v in arg_dict.items():
    if k not in opt_keys:
        continue
    if v is not None:
        opts += " --{} \"{}\" ".format(k.replace("_", "-"), v)

if arg_dict["wrap"] is not None:
    cmd = "sbatch {opts}".format(opts=opts)
else:
    cmd = "sbatch {opts} {extras}".format(opts=opts, extras=extras)

try:
    res = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    raise e

# Get jobid
res = res.stdout.decode()
try:
    m = re.search("Submitted batch job (\d+)", res)
    jobid = m.group(1)
    print(jobid)
except Exception as e:
    print(e)
    raise
