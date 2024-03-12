#!/bin/sh

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/DATA/pengchao/apps/anaconda/3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/DATA/pengchao/apps/anaconda/3/etc/profile.d/conda.sh" ]; then
        . "/DATA/pengchao/apps/anaconda/3/etc/profile.d/conda.sh"
    else
        export PATH="/DATA/pengchao/apps/anaconda/3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<


module load anaconda/anaconda.2020.02 && source activate deepkskit_0.1
ulimit -s unlimited

nohup python -u -m deepks iterate machines.yaml params.yaml systems.yaml scf_abacus.yaml >> log.iter 2> err.iter & 
