## Foodweb_thermal_asymmetries

**Data:**
1) Most data used in this paper can be found in the DATA folder. Some data (bird data from BirdLife International) was too large to report here, but can still be freely downloaded from BirdLife International as explained in the Methods section. 

**Data anlysis:**
1) Temp_data_endos.R:\
      –has all the code needed to prepare the data used for analysis as described in the Methods section. It requires the data        from the DATA folder.

2) Asymmetry_plots.R:\
      –has the code needed to reproduce Figure 2 of the main text. 

3) Ea_TL_analyses.R:\
      –has the code needed to reproduce Figure 4 of the main text. It requires data from the DATA folder.
      –the code 1) loads the data in Dell et al 2011, 2) finds all those ocurrences for which information on possible trophic       role exists (e.g., resource, grazer, omnivore), 3) extracts reported temperature sensitivities for those species, then         4) asks whether the temperature sensitivities increase or decrease with trophic level, by performing the analysis in two       different ways (as describved in the main text). The code is setup to reproduce the results in the main text. 

**Modeling:**
1) Model_Julia_Code:\
      –explains in pseudocode how the Julia code used to run the model in the main text works

2) Model_Julia_Code:\
      –has the code needed to reproduce Figure 3 of the main text and Figure S1 of the appendix\
      –to use, please download Julia from https://julialang.org/. Notice that all coding was done in Julia v1.0.0\
      –all Julia packages needed to run the code can be found in lines 8-12, and will also have to be downloaded to run the         script\
      –the code 1) defines the model in differential equations described in main text as a function, 2) then analyzes the           three scenarios described in the main text (only mortality rates are temperature-dependent, only attack rates are             temperature-dependent, and both parameters are simultaneously temperature-dependent). The code is already setup to             reproduce the results in the main text.
      –the code runs faster if pre-compiled. To do so, set the values in each for loop to i in 1:1 instead of 1:12/1:100 and         run once. Then run again with the provided values i in 1:12/1:100. Expected running time should be around 3-4 minutes on       a 36 core computer (Mac OS) for mortality and attack rates models, and 6-12 minutes for attack and mortality rates run.
