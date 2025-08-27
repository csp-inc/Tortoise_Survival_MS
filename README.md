# Tortoise_Survival_MS
This repository contains R code and data for replication of the analysis used in the Endangered Species Research manuscript 'An integrated model improves inferences about survival in the Mojave desert tortoise'

There are five files provided that are used to set up and run the JAGS models: 

'Tortoise_Survival_Functions.R' provides functions that write out the .txt JAGS code files and runs the analysis.

'Tortoise_Survival_Analysis.R' provides code to load the data and run the analysis.

'surv_data_telemetry_public.RDS', 'surv_data_markrecapture_public.RDS', and 'surv_data_integrated_public.RDS' are the datasets for each of the different models (known-fate, Cormack-Jolly-Seber, and integrated) formatted for the JAGS analysis. 

To protect sensitive tortoise locations, all site names have been removed and UTM locations have been randomly shifted from the true locations but kept in the same within-plot spatial relationship.

Running the CJS and Integrated models can take over 24 hours depending on computation resources.
