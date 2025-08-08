library(tidyverse)
library(argparser)
library(readxl)

# Create an argument parser
p = arg_parser("Prepare MVP lab files for GSR data tables")
p = add_argument(p, "--gsr-filelist", help = "Path to the tsv containing GSR file paths, chromosome, and md5sum")
p = add_argument(p, "--gsr-source-url", help = "URL of the original GSR file")
p = add_argument(p, "--metadata-file", help = "Path to the metadata file containing additional information")
p = add_argument(p, "--pha-mapping-file", help = "Path to the pha mapping file on dbGaP")
p = add_argument(p, "--contributor", help = "Email address of person running this script")
p = add_argument(p, "--output-dir", help = "Directory to save output files", default = "output")
argv = parse_args(p)

# Make sure the output directory does not exist.
if (dir.exists(argv$output_dir)) {
    stop("Output directory already exists: ", argv$output_dir)
}

# Read in GSR filelist.
gsr_files <- read_tsv(argv$gsr_filelist)
gsr_files %>% print(n=Inf)

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
    "gsr_source_url" = argv$gsr_source_url,
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

# Create the association files table.
gsr_files = gsr_files %>%
    rename(
        file_path=file,
        chromosome=chr,
    )
gsr_files

# Write out tables for loading into data tables.
dir.create(argv$output_dir)

# Write out the association analysis file.
association_analysis_file = file.path(argv$output_dir, "association_analysis.tsv")
write_tsv(tibble::enframe(association_analysis), association_analysis_file)

# Write out the files table.
gsr_files_file = file.path(argv$output_dir, "association_file.tsv")
write_tsv(gsr_files, gsr_files_file)
