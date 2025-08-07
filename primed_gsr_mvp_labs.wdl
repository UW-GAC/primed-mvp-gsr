version 1.0

workflow prep_gsr_mvp_lab {
    input {
        String gsr_file
        String metadata_file
        String contributor_email
        String pha_mapping_file="https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs002453/analyses/phs002453.v1.p1_GIA_pha_mapping.xlsx"
    }

    call prep_gsr_mvp_lab {
        input:
            gsr_file=gsr_file,
            metadata_file=metadata_file,
            pha_mapping_file=pha_mapping_file,
            contributor_email=contributor_email
    }

    meta {
          author: "Adrienne stilp"
          email: "amstilp@uw.edu"
    }
}

task prep_gsr_mvp_lab {
    input {
        String gsr_file
        String metadata_file
        String pha_mapping_file
        String contributor_email
    }

    command <<<
        set -e -o pipefail
        Rscript /usr/local/primed-mvp-gsr/prep_mvp_lab_files.R \
            --gsr-file ~{gsr_file} \
            --metadata-file ~{metadata_file} \
            --pha-mapping-file ~{pha_mapping_file} \
            --contributor ~{contributor_email} \
            --output-dir output
    >>>

    output {
        File association_analysis_file = "output/association_analysis.tsv"
        Array[File] association_files = glob("output/*.chr*.txt.gz")
    }

    runtime {
        docker: "uwgac/primed-mvp-gsr:0.0.1"
    }
}
