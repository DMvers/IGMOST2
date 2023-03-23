# Infant Gut Microbiota Over Space and Time V2

--INSTALLING--  
Tested on Ubuntu 18 and 21.

Utilised libraries on Ubuntu 18:
libsbml5-dev  
qt5-default  
libglpk-dev  

Utilised libraries on Ubuntu 21:
qtbase5-dev  
qtchooser  
qt5-qmake  
qtbase5-dev-tools  
libsbml5-dev  
libgsl-dev  
xml2  
libxml++2.6-dev  
libfreetype-dev  
libglpk-dev  

These can all be retrieved from standard repositories using, e.g., 'sudo apt-get install libsbml5-dev'

By default the software also utilises libpngwriter for drawing visualisations, which can be retrieved from https://github.com/pngwriter/pngwriter 
pngwriter can be utilised by this software by compiling pngwriter from source and placing the pngwriter.h file in the /include directory, renaming the libPNGwriter.a file to libpngwriter.a, and placing that in the /lib directory

These are suggestions - the code can be modified to, for example, use a different solver.

--DEFAULT RESOURCE FILES--  
Utilised bacterial models, all from AGORA version 1.03, originally retrieved from www.vmh.life:  
Bacteroides_vulgatus_ATCC_8482
Bifidobacterium_longum_infantis_ATCC_15697  
Bifidobacterium_longum_longum_ATCC_55813  
Clostridium_butyricum_DSM_10702
Collinsella_aerofaciens_ATCC_25986  
Eggerthella_sp_YY7918
Enterococcus_faecalis_OG1RF_ATCC_47077  
Escherichia_coli_SE11  
Eubacterium_hallii_DSM_3353
Gemella_morbillorum_M424
Haemophilus_parainfluenzae_T3T1
Lactobacillus_gasseri_ATCC_33323  
Parabacteroides_distasonis_ATCC_8503  
Propionibacterium_acnes_KPA171202
Roseburia_inulinivorans_DSM_16841
Rothia_mucilaginosa_DY_18
Ruminococcus_gnavus_ATCC_29149  
Staphylococcus_epidermidis_ATCC_12228
Streptococcus_oralis_Uo5
Veillonella_dispar_ATCC_17748  

The original versions that we used as a base for our changes can be found at https://github.com/DMvers/IGMOSTdatafiles

Newer versions may be available from www.vmh.life

Rename the model files to the following names for the default settings to work:
MODEL_Bvul.xml
MODEL_BlongInf.xml  
MODEL_BlongLong.xml  
MODEL_Cbut.xml
MODEL_Caer.xml
MODEL_EYY.xml
MODEL_Efae.xml
MODEL_Ecol.xml
MODEL_Ehal.xml
MODEL_Gmor.xml
MODEL_Hpar.xml
MODEL_Lgas.xml
MODEL_Pdis.xml
MODEL_Pacn.xml
MODEL_Rinu.xml
MODEL_Rmuc.xml
MODEL_Rgna.xml
MODEL_Sepi.xml
MODEL_Sora.xml
MODEL_Vdis.xml

Change the "name" characteristic on line 3 of the models to the following for the default settings to workBvul
Binf
Blong  
Cbut
Caer
EYY
Efae
Ecol
Ehal
Gmor
Hpar
Lgas
Pdis
Pacn
Rinu
Rmuc
Rgna
Sepi
Sora
Vdis


Now apply the changes described in S1 Table to these models to create the versions we used in our paper

--USE--  
Compile and run a simulation with the command "bash run_simulation.sh"
Please ensure there are no spaces in the directory name

The model can be run with different parameters by altering run_simulation.sh
Different bacterial models (any from www.vmh.life should be compatible) can be used by changing the references to xml files in simulation.cpp
