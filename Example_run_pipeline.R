# Practice:
library(targets)

# 1. Set the data, configuration, parameters in _targets.R
# It is recommended to use R project

# 2. Run the pipeline before rank-normalization
tar_make(names = c(data, dproc, dproc_csv))

# 3. After setting the configuration following the check results:
# Run the rest of the pipeline, target will skip the checking step, as it's already run
tar_make()

# 4. Check your temporarily saved data in targets hidden folder 
# "_targets/objects/", for example:
filtering_step<- tar_read(filtered_res)
# Now you can save the filtering results if needed
write.csv(filtering_step[['htn_incident']],'./Example_EN/Output/Filtering_htn.csv',row.names = F)

