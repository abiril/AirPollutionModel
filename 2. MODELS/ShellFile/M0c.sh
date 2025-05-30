#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=24:mem=128gb:ompthreads=24
#PBS -N M0c
#PBS -o /rds/general/user/air21/home/MODELS/MODEL/Outputs/
#PBS -e /rds/general/user/air21/home/MODELS/MODEL/Outputs/

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate RINLA-R4.3.3

cd $HOME

R CMD BATCH $HOME/MODELS/MODEL/Code/M0c.R $HOME/MODELS/MODEL/Outputs/M0c.ROut