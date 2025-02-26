## The OGSTM-BFM-Hg model

OGSTM-BFM-Hg is a transport-reaction model coupling mercury biogeochemistry in the ocean with the plankton-nutrients-detritus dynamics (BFM model) and hydrodynamic transport (OGSTM model). The model simulates the main marine Hg species (HgII, Hg0, MMHg, and DMHg) as well as other key processes of the Hg cycle such as bioaccumulation in phyto- and zooplankton, partitioning to detritus, sinking to the seabed, and outgassing to the atmosphere. The tool is written in Fortran 90 and presented in Rosati et al., (2022). It was developed under the Italian PRIN project ICCC (Impact of Climate Change on the biogeochemistry of Contaminants in the Mediterranean Sea) and improved as part of the NECCTON project. 

**References**

*Rosati, G., Canu, D., Lazzari, P., & Solidoro, C. (2022). Assessing the spatial and temporal variability of MeHg biogeochemistry and bioaccumulation in the Mediterranean Sea with a coupled 3D model. Biogeosciences, 19(February), 3663–3682. https://doi.org/doi.org/10.5194/bg-2022-14*


**Installation**

The shell script *downloader_ogstm_bfm.sh* automatically downloads the OGSTM and BFM codes from their repositories and switches to the branch containing the mercury dynamics (branch: neccton_WP8). Once the model directories have been downloaded, the *builder_ogstm_bfm.sh* can be executed to compile the model. These file are available at: https://github.com/inogs/ModelBuild/tree/neccton_WP8. 
For the procedure to be successful, a request to access the BFM repository must be made in advance to the BFM system team via this link https://docs.google.com/forms/d/e/1FAIpQLScI7N8AcvFxBeCD-EXwMXkQhgMwjhOLz3MYX8Kb47oPCXRv6w/viewform  

**Required data**

The model namelists, generated during compilation, must be copied in the run directory ‘~/wrkdir/MODEL’:

cp ${CODEDIR}/ogstm/ready_for_model_namelists/* .

The run directory must include a FORCINGS directory containing the physical fields for running the model, and RESTARTS and BC directories containing files for initial and boundary conditions of biogeochemical and mercury variables. A file of the domain, called “meshmask.nc”, is required to run the model. All files used are in the netCDF format.
