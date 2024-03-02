
The following numerical simulation is divided into 2 different parts.
The first is calculating the differnt datasets of the 2D thermodynamic model of 56 subduction zones used in the publication. 
The second part is combining all settings and calculated results into a self explainatory datastructur.

Model 1) Dataset 1 with Mg-sursassite
Model 2) Dataset 2 with Mg-sursassite
Model 6) Dataset 2 without Mg-sursassite

This numerical model is programmed to run on MacOS 10.14.6 and MATLAB R2021b.
The Project folder needs read and write permission for everyone to run Perple_X.
MacOS might require to execute /Final_Gies_et_al_2023/SZ_Model/vertex and /Final_Gies_et_al_2023/SZ_Model/werami before the first use with a rightclick>open


Instructions: 
1) Set your working and output directory in the run_SZ.m function
   line 9 base_dir='/change/to/your/path/Final_Gies_et_al_2023'; 

the total number of output files in the folder your/path/Final_Gies_et_al_2023/data_output_ext_int1km/SZ_results has to be 168.

2) Set the output direcory of run_SZ.m in the combine_results.m script.
   line 1 results_data=dir('/change/to/your/path/Final_Gies_et_al_2023/data_output_ext_int1km/SZ_results/');

3) Excecute in MATLAB run_SZ.m
this can be done in multiple MATLAB instances to reduce calculation time.

For paralell computation run the following line in multiple terminal windows.
make sure you change the all data paths "/change/to/your/" to your datapaths.

/change/to/your/MATLAB/path/Applications/MATLAB_R2021b.app/bin/matlab -batch "sznumin=(1:56);addpath(genpath('/change/to/your/path/Final_Gies_et_al_2023'));run_SZ(sznumin)" -sd "/change/to/your/path/Final_Gies_et_al_2023"

4) Excecute in MATLAB: addpath(genpath('/change/to/your/path/Final_Gies_et_al_2023')); 

5) Excecute in MATLAB combine_results.m to organize the model results in a datastructure.

Scripts for data extraction and data visualisation are available on request from Nils Gies
nils.gies@geo.unibe.ch
