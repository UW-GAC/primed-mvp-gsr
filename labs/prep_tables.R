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


# Prepare association analysis table.
association_analysis = c()

# Analysis identifiers.
association_analysis["gsr_source"] = "dbGaP"
association_analysis["gsr_source_url"] = argv$gsr_source_url
association_analysis["dbgap_analysis_accession"] = pha
association_analysis["pubmed_id"] = "39024449"
association_analysis["first_author"] = "Verma A"
association_analysis["publication_url"] = "https://pubmed.ncbi.nlm.nih.gov/39024449/"
association_analysis["consent_code"] = "NRES" # Is this correct? DUO term for unrestricted access.
association_analysis["upload_date"] = lubridate::today()
association_analysis["contributor_contact"] = argv$contributor # Pull from current AnVIL user?



# Extract trait information.
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
association_analysis["trait"] = trait_info$trait
association_analysis["trait_type"] = "quantitative"
association_analysis["trait_unit"] = trait_info$unit
association_analysis["trait_transformation"] = transformation
association_analysis["trait_definition"] = metadata$Description[metadata$Info == "Analyzed variable"]
association_analysis["covariates"] = metadata$Description[metadata$Info == "Covariates"]

# Variant information.
association_analysis["reference_assembly"] = "GRCh38"
association_analysis["n_variants"] = gsr_files$n_variants %>% as.integer() %>% sum(na.rm=T)
association_analysis["genotyping_platform"] = "Affymetrix MVP 1.0 Genotyping Array"
association_analysis["genotyping_technology"] = "genome-wide array"
association_analysis["is_imputed"] = TRUE
association_analysis["imputation_reference_panel"] = "Other"
association_analysis["imputation_reference_panel_detail"] = "African Genome Resources panel and 1000 Genomes Project (1000G) Phase 3 Version 5"
association_analysis["imputation_quality_filter"] = 0.3
association_analysis["cohorts"] = "MVP"

# Sample information.
n_samp = metadata %>%
    filter(Info == "Sample size") %>%
    mutate(
        tmp = str_match(Description, "Total Sample Size=(.+?)$")[, 2],
    ) %>%
    pull(tmp) %>%
    as.integer()
population_label = metadata %>%
    filter(Info == "Sample population") %>%
    mutate(
        tmp = str_match(Description, "^(.+?) \\(GIA\\)$")[, 2],
    ) %>%
    pull(tmp)
association_analysis["n_samp"] = n_samp
association_analysis["n_effective"] = n_samp
association_analysis["population_descriptor"] = "Genetically inferred ancestry"
association_analysis["population_labels"] = population_label
association_analysis["is_meta_analysis"] = str_detect(argv$gsr_source_url, ".META.")

# Analysis information.
association_analysis["analysis_method"] = "LMM"
association_analysis["analysis_software"] = "SAIGE"

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
association_analysis = tibble::enframe(association_analysis, name="field", value="value")
write_tsv(association_analysis, association_analysis_file)

# Write out the files table.
gsr_files_file = file.path(argv$output_dir, "association_file.tsv")
write_tsv(gsr_files, gsr_files_file)
