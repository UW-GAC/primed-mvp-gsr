library(tools)
library(tidyverse)
library(glue)
library(argparser)

# Create an argument parser
p = arg_parser("Prepare MVP lab files for GSR data tables")
p = add_argument(p, "--gsr-file", help = "Path to the input file containing GSR data")
p = add_argument(p, "--output-directory", help = "Directory (or cloud path) where output files should be saved")
argv = parse_args(p)

print(argv)

# Check if the output directory already exists.
output_directory = argv$output_directory
if (!str_starts(output_directory, "gs://") & dir.exists(output_directory)) {
    stop("Output directory already exists: ", output_directory)
}

# Download the file.
gsr_file = basename(argv$gsr_file)
download.file(argv$gsr_file, destfile=gsr_file)

# Read in GSR file.
gsr_data <- read_tsv(gsr_file)
gsr_data %>% print(n=6)

# Rename columns as necessary.
#  [1] "SNP_ID"      "chrom"       "pos"         "ref"         "alt"
#  [6] "ea"          "af"          "num_samples" "beta"        "sebeta"
# [11] "pval"        "r2"          "q_pval"      "i2"          "direction"
gsr_data = gsr_data %>%
    rename(
        rsID=SNP_ID,
        chromosome=chrom,
        position=pos,
        effect_allele=ea,
        # Note that ref and alt are not dbSNP ref/alt, but a different standard.
        # See eLetters here: https://www.science.org/doi/10.1126/science.adj1182
        other_allele=ref,
        effect_allele_freq=af,
        p_value=pval,
        se=sebeta,
        n_samp=num_samples,
        imputation_quality_score=r2,
        heterogeneity_p_value=q_pval,
        heterogeneity_I2=i2,
    ) %>%
    mutate(
        # Spot checks imply that it is the forward strand.
        strand="forward"
    ) %>%
    select(
        -alt
    )

tempdir = tempfile()
dir.create(tempdir)
# Get basename of the files so we can append chromosome.
file_base = gsr_file %>%
    basename() %>%
    str_replace(".txt.gz", "")
# Get the final directory and remove the trailing slash if it exists.
output_directory = argv$output_directory
output_directory = str_replace(output_directory, "/$", "")

# Split the GSR data by chromosome and store relevant information.
gsr_data = gsr_data %>%
    mutate(chr = chromosome) %>%
    group_by(chr) %>%
    nest() %>%
    mutate(
        # Number of variants.
        n_variants = map_int(data, ~ nrow(.x)),
        basename = paste0(file_base, ".chr", chr, ".txt.gz"),
        tempfile = file.path(tempdir, basename),
        file = file.path(output_directory, basename),
        file_type = "data",
    )

# Write out the temporary files and get md5sums.
for (i in seq_along(gsr_data$chr)) {
    write_tsv(gsr_data$data[[i]], gsr_data$tempfile[i], na = "", append = FALSE)
}

gsr_files = gsr_data %>%
    select(-data)

# Create a data dictionary.
data_dictionary = c(
    rsID = "rs identifier",
    chromosome = "the chromosome that the variant is located on",
    position = "the base pair location of the variant",
    other_allele = "other allele of the variant",
    effect_allele = "effect allele of the variant",
    effect_allele_freq = "frequency of the effect allele",
    n_samp = "number of samples for this variant",
    beta = "estimated effect size",
    se = "standard error ofbetra",
    p_value = "p-value",
    # Pulled description from GWAMA documentation.
    # https://genomics.ut.ee/en/tools
    imputation_quality_score = "Imputation r2 for the variant",
    heterogeneity_p_value = "For meta-analysis, Cochran's heterogeneity statistic's p-value",
    heterogeneity_I2 = "For meta-analysis, heterogeneity index I2 by Higgins et al 2003",
    direction = "For meta-analysis, the direction of the effect size in each group",
    strand = "DNA strand designation"
)
data_dictionary = data_dictionary[names(gsr_data$data[[1]])]
if (any(is.na(data_dictionary))) {
    stop("Data dictionary is missing some entries: ", paste(names(gsr_files)[is.na(data_dictionary)], collapse=", "))
}
data_dictionary = enframe(data_dictionary, name="field", value="value")

stopifnot(all(names(gsr_data$data[[1]]) == data_dictionary$field))

# Write out the data dictionary and save it
dd_file = file.path(tempdir, paste0(file_base, "_DD.tsv"))
write_tsv(data_dictionary, dd_file)
gsr_files = gsr_files %>%
    mutate(
        chr=as.character(chr),
        n_variants=as.character(n_variants),
    ) %>%
    bind_rows(tibble(
        tempfile = dd_file,
        chr = "None",
        file = file.path(output_directory, paste0(file_base, "_DD.tsv")),
        file_type = "data dictionary",
    ))

# Calculate md5sums.
gsr_files = gsr_files %>%
    mutate(md5sum = unname(md5sum(tempfile)))

# Save the file in the final location.
# If cloud path is provided, copy files to cloud storage.
if (str_starts(output_directory, "gs://")) {
    "gsutil -m cp {tempdir}/* {output_directory}" %>%
        glue() %>%
        system(intern = TRUE)

} else {
    # If not a cloud path, just move the files.
    file.rename(tempdir, output_directory)
}

# Write the list of files out.
gsr_files = gsr_files %>%
    select(file, chr, n_variants, md5sum, file_type)

print(gsr_files)
write_tsv(gsr_files, "gsr_files.tsv")
