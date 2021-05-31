To store code and datasets for Chapter 6: Multiple Imputation for Multivariate Missing Data: the Joint Modeling Approach

Example6.2.R uses the dataset colorado_ami_chf_pb.txt (the hosptial performance data). Example6.2_t_imputation.R applied the trivariate-t model for imputation. 

Example6.3.R uses the dataset paper1_dataset.xls (from Chapter 4).

Example6.4.R uses the dataset sbirth.txt.

Example6.7.R is a simulation program. It calls the WinBUGS program logitmodel_missing_bivariate_forR.txt.

Example6.8.R is a simulation program. It calls the WinBUGS program logitmodel_missinginteraction_forR.txt.

Example6.9.R uses the dataset week30.txt. It calls the WinBUGS program mixturemodel_complete.txt and mixturemodel_impute.txt.

For R programs that call WinBUGS (Examples 6.7.R, Example6.8.R, and Example6.9.R), the software WinBUGS is required to be installed. The simulation usually runs a long time. In the future we will provide the R NIMBLE version which does not require the WinBUGS software.
