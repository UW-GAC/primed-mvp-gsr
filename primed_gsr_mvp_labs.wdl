version 1.0

import "https://raw.githubusercontent.com/UW-GAC/primed-file-checks/refs/heads/main/validate_gsr_model.wdl" as validate_gsr_model

workflow prep_gsr_mvp_lab {
    input {
        String gsr_file
        String metadata_file
        String contributor_email
        String output_directory
        String pha_mapping_file="https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs002453/analyses/phs002453.v1.p1_GIA_pha_mapping.xlsx"
        String model_url
        String workspace_name
        String workspace_namespace
        Boolean overwrite = false
        Boolean import_tables = false
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

    call validate_gsr_model.validate_gsr_model {
        input:
            table_files = prep_gsr_tables.tables_file,
            model_url = model_url,
            workspace_name = workspace_name,
            workspace_namespace = workspace_namespace,
            overwrite = overwrite,
            import_tables = overwrite
    }

    output {
        File validation_report = validate_gsr_model.validation_report
        Array[File] tables = validate_gsr_model.tables
        String? md5_check_summary = validate_gsr_model.md5_check_summary
        File? md5_check_details = validate_gsr_model.md5_check_details
        String? data_report_summary = validate_gsr_model.data_report_summary
        File? data_report_details = validate_gsr_model.data_report_details
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

        cat '{"association_analysis": "","association_file": ""}' %>% "tables.txt"
    >>>

    output {
        File association_analysis_file = "output/association_analysis.tsv"
        File association_file = "output/association_file.tsv"
        Map[String, File] tables_file = read_map("tables.txt")

    }

    runtime {
        docker: "uwgac/primed-mvp-gsr:0.0.1"
    }

}
