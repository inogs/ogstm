This is an ensuite to run in-water RTM

1. Copy OASIM files from SCRATCH! Check if that's also in the archive:

	/gpfs/scratch/userexternal/plazzari/OASIM_HF_15m/DATA

2. Use unzip.sh along with create_lista_date.py to extract the *_0m.nc files you need for RTM - it takes time!

3. Modify the create_ave.py in order to save data for every 15 minutes, for each wl and variable (Ed, Es) separately
   You use it to restructure the input for the profiler.py

4. use create_ave.sh with lista_date.txt (from 01012012 to 31122017) and with create_ave.py

5. export MASKFILE=$PWD/meshmask.nc #/galileo/home/userexternal/eterzic0/CODE/ogstm/preproc/RADIATIVE_INWATER/meshmask.nc 

6. export PYTHONPATH=/galileo/home/userexternal/eterzic0/BIT_SEA_mod/bit.sea/:$PYTHONPATH

7. run the test.sh


7. run the profiler.py (change the path of VardescriptorB.xml accordingly)

8. run job.spawn.slurm

9. do the postproc


#source ../../compilers/machine_modules/galileo.intel 

#make -f Makefile_IOP compute_IOP.xx