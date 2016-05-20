# init parallel option
parallel = TRUE

# Process the data
source("./scripts/data_processing/reformat_salix.r")
source("./scripts/data_processing/split_data.r")

# Get the results
source("./scripts/results/1_metaweb_sampling.r")
source("./scripts/results/2_example_pair.r")
source("./scripts/results/3_map_pair.r")
source("./scripts/results/4_metaweb_holes.r")
source("./scripts/results/5_map_all.r")
source("./scripts/results/table_all_models.r")