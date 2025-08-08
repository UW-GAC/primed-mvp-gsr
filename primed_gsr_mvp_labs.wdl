version 1.0

workflow prep_gsr_mvp_lab {
    input {
        String gsr_file
        String metadata_file
        String contributor_email
        String output_directory
        String pha_mapping_file="https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs002453/analyses/phs002453.v1.p1_GIA_pha_mapping.xlsx"
        Int? mem_gb
    }

    call prep_gsr_files {
        input:
            gsr_file=gsr_file,
            output_directory=output_directory,
            mem_gb=mem_gb
    }

    call prep_gsr_tables {
        input:
            gsr_filelist=prep_gsr_files.gsr_files,
            gsr_source_url=gsr_file,
            metadata_file=metadata_file,
            pha_mapping_file=pha_mapping_file,
            contributor_email=contributor_email
    }

    output {
        File association_analysis_file = prep_gsr_tables.association_analysis_file
        File association_files_file = prep_gsr_tables.association_file
    }

    meta {
          author: "Adrienne stilp"
          email: "amstilp@uw.edu"
    }
}

task prep_gsr_files {
    input {
        String gsr_file
        String output_directory
        Int mem_gb = 16
    }

    command <<<
        set -e -o pipefail
        Rscript /usr/local/primed-mvp-gsr/labs/prep_gsr_files.R \
            --gsr-file ~{gsr_file} \
            --output-directory ~{output_directory}
    >>>

    output {
        File gsr_files = "gsr_files.tsv"
    }

    runtime {
        docker: "uwgac/primed-mvp-gsr:0.0.1"
        memory: "~{mem_gb} GB"
    }
}


task prep_gsr_tables {
    input {
        File gsr_filelist
        String metadata_file
        String pha_mapping_file
        String contributor_email
        String gsr_source_url
    }

    command <<<
        set -e -o pipefail
        Rscript /usr/local/primed-mvp-gsr/labs/prep_tables.R \
            --gsr-filelist ~{gsr_filelist} \
            --gsr-source-url ~{gsr_source_url} \
            --metadata-file ~{metadata_file} \
            --pha-mapping-file ~{pha_mapping_file} \
            --contributor ~{contributor_email}
    >>>

    output {
        File association_analysis_file = "output/association_analysis.tsv"
        File association_file = "output/association_file.tsv"
    }

    runtime {
        docker: "uwgac/primed-mvp-gsr:0.0.1"
    }

}
