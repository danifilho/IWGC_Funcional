configfile: "config.yaml"
import os

configfile: "config.yaml"

species_name = config['species_name']
volume_name = config['volume_name']

# Extract database information
mmseqs_databases = config['mmseqs_databases']

# Create lists of database names and paths
database_names = [db['name'] for db in mmseqs_databases]
database_paths = [db['path'] for db in mmseqs_databases]

# Create a dictionary mapping database names to paths
database_dict = {db['name']: db['path'] for db in mmseqs_databases}

# Identify the NCBI database
ncbi_database = next((db for db in mmseqs_databases if db['name'] == 'ncbi'), None)
if ncbi_database is None:
    raise ValueError("NCBI database not found in config.yaml")
ncbi_database_name = ncbi_database['name']
ncbi_database_path = ncbi_database['path']
ncbi_data_file = ncbi_database['ncbidata']  # The .prot file

# Optional: Print statements for debugging
print("Species Name:", species_name)
print("Database Names:", database_names)
print("Database Paths:", database_paths)
print("Database Dict:", database_dict)
print("NCBI Database Name:", ncbi_database_name)
print("NCBI Database Path:", ncbi_database_path)
print("NCBI Data File:", ncbi_data_file)

def match_database(clean_database_name):
    try:
        index = database_names.index(clean_database_name)
        database = database_paths[index]
        print(f"Matched database {clean_database_name} to path {database}")
        return database
    except ValueError as e:
        print(f"Error: Database name {clean_database_name} not found in paths {database_paths}.")
        raise e

def gather_inputs(species_name):
    mmseqs_outputs = expand("mmseqs_output/{species_name}_{database_name}_best_hit.out",
                            species_name=species_name, database_name=database_names)
    other_outputs = [
        f"hmmer_outputs/{species_name}.domtblout",
        f"multiloc_outputs/{species_name}_ml2.txt",
        f"sigtarp_outputs/{species_name}_summary.targetp2",
        f"sigtarp_outputs/{species_name}.gff3",
        f"iprscan_output/{species_name}.iprscan.tsv",
        "databases/current/GenomesDatabase.v1.prot"
    ]
    return mmseqs_outputs + other_outputs

rule all:
    input:
        f"functional_outputs/functional_annotation.tsv"

rule agat:
    input:
        final_filtering_gff="inputs/{config['final_filtering_gff']}"
    output:
        agat_output="agat_outputs/{species_name}.gff"
    singularity:
        "images/agat.sif"
    shell:
        "agat_sp_keep_longest_isoform.pl -gff {input} -o {output}"

rule gff_read:
    input:
        agat_output="agat_outputs/{species_name}.gff",
        fasta_file="inputs/{config['fasta_file']}"
    output:
        aa="gff_read_output/{species_name}.aa",
        cds="gff_read_output/{species_name}.cds"
    params:
        volume_name=config['volume_name']
    singularity:
        "images/gffread.sif"
    shell:
        """
        gffread {params.volume_name}/{input.agat_output} -g {params.volume_name}/{input.fasta_file} -J -S -y {params.volume_name}/{output.aa}
        gffread {params.volume_name}/{input.agat_output} -g {params.volume_name}/{input.fasta_file} -J -S -x {params.volume_name}/{output.cds}
        """

# Rule to run MMseqs2 for each database
rule mmseqs2:
    input:
        aa="gff_read_output/{species_name}.aa",
        database=lambda wildcards: database_dict[wildcards.database_name]
    output:
        ofmt6="mmseqs_output/{species_name}_{database_name}.ofmt6"
    params:
        volume_name=volume_name
    shell:
        """
        echo "Running mmseqs2 with input: {input.aa} and database: {input.database}"
        singularity exec -B {params.volume_name}:/data {params.volume_name}/images/mmseqs2.sif /bin/bash -c "cd /data && mmseqs easy-search {input.aa} {input.database} {output.ofmt6} /tmp --threads 20 --format-mode 0 --start-sens 2 -s 7 --sens-steps 3"
        """

# Rule to get best hits from MMseqs outputs
rule best_hits:
    input:
        ofmt6="mmseqs_output/{species_name}_{database_name}.ofmt6"
    output:
        best_hits_output="mmseqs_output/{species_name}_{database_name}_best_hit.out"
    params:
        volume_name=volume_name
    shell:
        """
        sort -k12 -t $'\t' -nr {input.ofmt6} | awk -F "\t" ' ! a[$1]++ && $11 < 1e-5' > {output.best_hits_output}
        """

rule interproscan:
    input:
        #gff_output=rules.best_hits.output,
        aa="gff_read_output/{species_name}.aa"
    output:
        iprscan_output="iprscan_output/{species_name}.iprscan.tsv"
    params:
        volume_name=config['volume_name']
    shell:
        """
        singularity exec -B {params.volume_name}/databases/interproscan-5.67-99.0/data:/opt/interproscan/data -B {params.volume_name}:/data {params.volume_name}/images/iprscan.sif /opt/interproscan/interproscan.sh --formats GFF3 TSV  --goterms --pathways --iprlookup --input /data/{input.aa} --cpu 30 --output-file-base /data/iprscan_output/{species_name}.iprscan
        """

rule prot_list:
    input:
        inteproscan_output=rules.interproscan.output,
        aa="gff_read_output/{species_name}.aa",
        iprscan_output="iprscan_output/{species_name}.iprscan.tsv"
    output:
        multigo="multiloc_outputs/{species_name}_multilocGo.out",
        prot_list="multiloc_outputs/{species_name}_prot.list"
    params:
        volume_name=config['volume_name']
    shell:
        """
        awk -F"\t" '$14 ~ /GO/ {{print $1,$14}}' {params.volume_name}/{input.iprscan_output} | sed 's/|/ /g' > {output.multigo}
        egrep \> {params.volume_name}/{input.aa} | cut -b 2- | awk '{{print $1}}' > {output.prot_list}
        mkdir splits
        cd splits
        split -n l/20 {params.volume_name}/{output.prot_list} #20 is the number of threads
        """

rule multiloc:
    input:
        prot_list_output=rules.prot_list.output,
        multiloc_script="scripts/{config['multiloc_script']}",
        aa="gff_read_output/{species_name}.aa",
        multigo="multiloc_outputs/{species_name}_multilocGo.out"
    output:
        multiloc_output="multiloc_outputs/{species_name}_ml2.txt"
    params:
        volume_name=config['volume_name']
    shell:
        """
        cd splits
        declare -a pid_array
        for x in x{{a..z}}{{a..z}}; do
            mkdir $x"_dir";
            {params.volume_name}/{input.multiloc_script} $x {input.aa} {input.multigo} $x"_dir" {params.volume_name} & pid_array+=("$!");
        done
            for pid in "${{pid_array[@]}}"; do
            wait $pid
        done
        cat x*/ml2.out > {params.volume_name}/{output.multiloc_output}
        cd -
        rm -rf splits
        """


rule hmmer:
    input:
        #gff_read_output=rules.multiloc.output,
        aa=f"gff_read_output/{species_name}.aa"
    output:
        hmmer_output="hmmer_outputs/{species_name}.domtblout"
    params:
        volume_name=config['volume_name']
    singularity:
        "images/hmmer3.sif"
    shell:
        """
        hmmscan --cpu 4 --domtblout {params.volume_name}/{output.hmmer_output} {params.volume_name}/databases/Pfam-A.hmm {params.volume_name}/{input.aa}
        """

rule sigtarp:
    input:
        rules.hmmer.output,
        aa="gff_read_output/{species_name}.aa"
    output:
        signalp_output="sigtarp_outputs/{species_name}.gff3",
        targetp_output="sigtarp_outputs/{species_name}_summary.targetp2"
    params:
        volume_name=config['volume_name']
    singularity:
        "images/sigtarp.sif"
    shell:
        """
        cd sigtarp_outputs
        signalp -fasta {params.volume_name}/{input.aa} -org euk -format short -gff3
        targetp -fasta {params.volume_name}/{input.aa} -org pl -format short
        cd -
        """

mmseqs_best_hits = [f"mmseqs_output/{species_name}_{db_name}_best_hit.out" for db_name in database_names]

rule generate_tsv:
    input:
        sigtar_signalp=f"sigtarp_outputs/{species_name}.gff3",
        sigtar_targetp=f"sigtarp_outputs/{species_name}_summary.targetp2",
        mmseqs=mmseqs_best_hits,
        iprs=f"iprscan_output/{species_name}.iprscan.tsv",
        hmmer=f"hmmer_outputs/{species_name}.domtblout",
        multiloc=f"multiloc_outputs/{species_name}_ml2.txt",
        ncbi_database=ncbi_data_file
    output:
        output_file=f"functional_outputs/functional_annotation.tsv"
    params:
        volume_name=volume_name,
        outputs_list=f"functional_outputs/{species_name}_outputs_list"
    run:
        import os
        # Create the outputs_list file dynamically
        with open(params.outputs_list, 'w') as f:
            for item in input.mmseqs:
                f.write(f"{params.volume_name}/{item}\n")
            f.write(f"{params.volume_name}/{input.iprs}\n")
            f.write(f"{params.volume_name}/{input.hmmer}\n")
            f.write(f"{params.volume_name}/{input.multiloc}\n")
            f.write(f"{params.volume_name}/{input.sigtar_signalp}\n")
            f.write(f"{params.volume_name}/{input.sigtar_targetp}\n")
            f.write(f"{params.volume_name}/{input.ncbi_database}\n")
        shell(f"""
            python {params.volume_name}/scripts/functional_merge.py --filelist {params.volume_name}/{params.outputs_list} --output_dir {params.volume_name}/functional_outputs --species_name {species_name}
            """)
        # Optionally, remove the outputs_list file after use
#        os.remove(params.outputs_list)

