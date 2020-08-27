# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/miniconda3-4.5.4-hivczbzklvoccmuifapprxz7humnmn5c/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/miniconda3-4.5.4-hivczbzklvoccmuifapprxz7humnmn5c/etc/profile.d/conda.sh" ]; then
        . "/usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/miniconda3-4.5.4-hivczbzklvoccmuifapprxz7humnmn5c/etc/profile.d/conda.sh"
    else
        export PATH="/usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/miniconda3-4.5.4-hivczbzklvoccmuifapprxz7humnmn5c/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

export PATH=/home/xwy21/src/MACS-1.4.2/bin:$PATH
export PATH=/home/xwy21/src/cellranger-atac-1.2.0:$PATH
export PATH=/home/xwy21/src/STAR/source:$PATH
export PATH=/home/xwy21/src/subread-2.0.1-Linux-x86_64/bin:$PATH
export PATH=/home/xwy21/src/genometools-1.6.1/bin:$PATH
export PATH=/home/xwy21/.local/bin:$PATH

module load libxml2-2.9.8-gcc-5.4.0-sy2r4k7
module load r-4.0.2-gcc-5.4.0-xyx46xb
