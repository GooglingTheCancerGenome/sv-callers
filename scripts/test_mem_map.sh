#!/bin/bash
#SBATCH --time-min=1
#SBATCH --mem=100
#sysctl -w vm.max_map_count=131072 &&
sysctl -a |grep vm.max_map 2>&1

