#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=60:mem=256gb:ompthreads=60
#PBS -N M1b
#PBS -o /rds/general/user/air21/home/MODELS/MODEL/Outputs/
#PBS -e /rds/general/user/air21/home/MODELS/MODEL/Outputs/

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate RINLA-R4.3.3

cd $HOME

R CMD BATCH $HOME/MODELS/MODEL/Code/M1b.R $HOME/MODELS/MODEL/Outputs/M1b.ROut