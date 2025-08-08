library(tidyverse)
library(glue)
library(argparser)

# Create an argument parser
p = arg_parser("Prepare MVP lab files for GSR data tables")
p = add_argument(p, "--gsr-file", help = "Path to the input file containing GSR data")
p = add_argument(p, "--cloud-path", help = "Location of output files in cloud storage")
p = add_argument(p, "--output-directory", help = "Directory to save output files", default = "output")

argv = parse_args(p)
print(argv)

# Check if the output directory already exists.
output_directory = argv$output_directory
if (dir.exists(output_directory)) {
    stop("output directory exists")
}
dir.create(output_directory)

# Download the file.
gsr_file = basename(argv$gsr_file)
download.file(argv$gsr_file, destfile=gsr_file)

# Read in GSR file.
gsr_data <- read_tsv(gsr_file)
gsr_data %>% print(n=6)

# Split the file by chromosome
gsr_data_split = gsr_data %>%
    group_by(chrom) %>%
    group_split()
lapply(gsr_data_split, print, n=6)


# Write out association files.
file_pattern = file.path(output_directory, basename(gsr_file)) %>%
    str_replace(".txt.gz", ".chr{chrom}.txt.gz")

for (x in gsr_data_split) {
    chrom = unique(x$chrom)
    print(chrom)
    stopifnot(length(chrom) == 1)
    outfile = file_pattern %>% glue(chrom=chrom)
    print(outfile)
    write_tsv(x, outfile, na = "", append = FALSE)
}

# If cloud path is provided, copy files to cloud storage.
cloud_path = argv$cloud_path
if (!is.na(cloud_path)) {
    # Remove trailing slash if it exists.
    cloud_path = str_replace(cloud_path, "/$", "")
    "gsutil -m cp {output_directory}/* {cloud_path}" %>%
        glue() %>%
        system(intern = TRUE)
    final_files = file.path(cloud_path, list.files(output_directory))
    unlink(output_directory, recursive = TRUE)
} else {
    final_files = file.path(output_directory, list.files(output_directory))
}

# Save the final files to a text file.
x = tibble(file = final_files) %>%
    # Extract the chromosome from the file name using a regular expression.
    mutate(
        chrom = str_match(file, "chr(.+?).txt.gz")[, 2]
    )
print(x)
write_tsv(x, "gsr_files.txt")
