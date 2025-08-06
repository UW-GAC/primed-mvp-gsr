library(tidyverse)
library(argparser)
library(readxl)

# Create an argument parser
p = arg_parser("Prepare MVP lab files for GSR data tables")
# Either pha or gsr_file and metadata_file are required.
# p = add_argument(p, "--pha", help = "PHA for this GSR analysis", type="integer")
p = add_argument(p, "--gsr-file", help = "Path to the input file containing GSR data")
p = add_argument(p, "--metadata-file", help = "Path to the metadata file containing additional information")
# Other info.
p = add_argument(p, "--output-dir", help = "Directory to save output files", default = "output")
p = add_argument(p, "--pha-mapping-file", help = "Path to the pha mapping file on dbGaP")
p = add_argument(p, "--contributor", help = "Email address of person running this script")
argv = parse_args(p)

# # For testing.
# argv = list(
#     gsr_file="tmp/tmp_subset.txt.gz",
#     metadata_file="tmp/tmp.metadata.txt",
#     pha_mapping_file="tmp/pha_mapping.xlsx",
#     output_dir="output"
# )
# print(argv)

# Note: reading directly from the url of the file only returns a subset of the data.
# Possibly related to this GitHub issue:
# https://github.com/tidyverse/readr/issues/1555
# Instead we will have to download the file first, and then read it.
is_remote = (
    str_starts(argv$gsr_file, "https://") |
    str_starts(argv$gsr_file, "http://") |
    str_starts(argv$gsr_file, "ftp://")
)
print(is_remote)
if (is_remote) {
    gsr_source_url = argv$gsr_file
    gsr_file = tempfile(fileext = ".txt.gz")
    download.file(argv$gsr_file, destfile=gsr_file)
} else {
    gsr_source_url = NA
    gsr_file = argv$gsr_file
}

# Read in GSR file.
gsr_data <- read_tsv(gsr_file)
gsr_data %>% print(n=6)

# Read in metadata file.
metadata <- read_tsv(argv$metadata_file)
metadata %>% print(n=Inf)

# Read in the PHA mapping file.
is_remote = (
    str_starts(argv$pha_mapping_file, "https://") |
    str_starts(argv$pha_mapping_file, "http://") |
    str_starts(argv$pha_mapping_file, "ftp://")
)
if (is_remote) {
    tmpfile = tempfile(fileext = ".txt")
    download.file(argv$pha_mapping_file, destfile=tmpfile)
    pha_mapping <- read_excel(tmpfile)
} else {
    pha_mapping <- read_excel(argv$pha_mapping_file)
}
names(pha_mapping) = c("pha", "analysis_name")

# Get the pha for this GSR.
pha = pha_mapping %>%
    filter(analysis_name == metadata$Description[metadata$Info == "Title of analysis"]) %>%
    pull(pha)
print(pha)
stopifnot(length(pha) == 1)

# Get the transformation
analyzed_variable = metadata %>%
    filter(Info == "Analyzed variable") %>%
    pull(Description)
print(analyzed_variable)
stopifnot(length(analyzed_variable) == 1)

# Extract the part in parentheses; these have more information about the units and transformations.
trait_info = metadata %>%
    filter(Info == "Analyzed variable") %>%
    mutate(
        tmp1 = str_match(Description, "(.+?) \\((.+)\\)")[, 2],
        tmp2 = str_match(Description, "(.+?) \\((.+)\\)")[, 3]
    ) %>%
    separate_wider_delim(tmp2, names=c("abbreviation", "unit", "dist", "transformation"), delim=", ", too_few="align_end") %>%
    mutate(
        trait = paste(dist, tmp1)
    )
stopifnot(nrow(trait_info) == 1)

if ("inv-norm transformed" == trait_info$transformation) {
    transformation = "inverse normal"
} else {
    stop("Unknown transformation in analyzed variable: ", analyzed_variable)
}

# Prepare association analysis table.
association_analysis = c(
    # General info.
    "gsr_source" = "dbGaP",
    "gsr_source_url" = gsr_source_url,
    "dbgap_analysis_accession" = pha,
    "pubmed_id" = "39024449",
    "first_author" = "Verma A",
    "publication_url" = "https://pubmed.ncbi.nlm.nih.gov/39024449/",
    "consent_code" = "NRES", # Is this correct? DUO term for unrestricted access.
    "upload_date" = lubridate::today(),
    "contributor_contact" = argv$contributor_contact, # Pull from current AnVIL user?
    # Trait info
    trait = trait_info$trait,
    trait_type = "quantitative",
    trait_unit = trait_info$unit,
    trait_transformation = transformation,
    trait_definition = metadata$Description[metadata$Info == "Analyzed variable"],
    covariates = metadata$Description[metadata$Info == "Covariates"]
)
print(association_analysis)

# Split by chromosome
gsr_data_split = gsr_data %>%
    group_by(chrom) %>%
    group_split()
lapply(gsr_data_split, print, n=6)

# Create the association files table.
# TODO.

# Write out association files.
dir.create(argv$output_dir)
for (i in seq_along(gsr_data_split)) {
    chrom_data = gsr_data_split[[i]]
    chrom = unique(chrom_data$chrom)
    print(chrom)
    stopifnot(length(chrom) == 1)
    outfile = file.path(argv$output_dir, basename(argv$gsr_file)) %>%
        str_replace(".txt.gz", paste0(".chr", chrom, ".txt.gz"))
    print(outfile)
    write_tsv(chrom_data, outfile, na = "NA", append = FALSE)
}

# Write out the association analysis file.
association_analysis_file = file.path(argv$output_dir, "association_analysis.tsv")
write_tsv(tibble::enframe(association_analysis), association_analysis_file)
