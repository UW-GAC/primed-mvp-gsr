version 1.0

workflow prep_gsr_mvp_lab {
    input {
        String gsr_file
        String metadata_file
        String contributor_email
        String pha_mapping_file="https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs002453/analyses/phs002453.v1.p1_GIA_pha_mapping.xlsx"
        Int? mem_gb
    }

    call prep_gsr_mvp_lab {
        input:
            gsr_file=gsr_file,
            metadata_file=metadata_file,
            pha_mapping_file=pha_mapping_file,
            contributor_email=contributor_email,
            mem_gb=mem_gb
    }

    call prep_gsr_mvp_tables {
        input:
            association_file_paths=prep_gsr_mvp_lab.association_files
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
        Int mem_gb = 16
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
        memory: "~{mem_gb} GB"
    }
}

task prep_gsr_mvp_tables {
    input {
        Array[String] association_file_paths
    }

    command <<<
        echo ~{sep=" " association_file_paths}
    >>>

    runtime {
        docker: "uwgac/primed-mvp-gsr:0.0.1"
    }
}
