#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=64:mem=256gb:ompthreads=64
#PBS -N M4b_Laplace
#PBS -o /rds/general/user/air21/home/MODELS/MODEL/Outputs/
#PBS -e /rds/general/user/air21/home/MODELS/MODEL/Outputs/

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate RINLA-R4.3.3

cd $HOME

R CMD BATCH $HOME/MODELS/MODEL/Code/M4b_Laplace.R $HOME/MODELS/CHOSEN_ONE/Outputs/M4b_Laplace.ROut