#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.immunoinformatics_allele_table_path = '/home/pathinformatics/immunoinformatics_platform/host/host-data/allele-frequency-processed-tables/all_allele_frequencies.tsv'
params.tcrpmhc_template_file = '/home/pathinformatics/immunoinformatics_platform/host/host-data/tcrpmhc_templates/TCRpMHC_template_0.fasta'

/*params.protein_file = '/home/pathinformatics/example_data_for_nextflow/fasta_files/proteins/test.fasta'*/
params.protein_file = '/home/pathinformatics/example_data_for_nextflow/fasta_files/proteins/alphavirus_protein_multiseq.fasta'
params.epitope_output_folder = '/home/pathinformatics/epitope_outputs'

params.cdhit_similarity_threshold = 0.95
params.alphafold_pdb_folder = '/home/pathinformatics/epitope_outputs/alphafold_predictions/*/*.pdb'

params.netmhcpan_chunk_size = 500

params.cdhit_input_proteins = "no"
params.bepipred = "no"
params.epidope = "no"
params.netmhcpani = "no"
params.netmhcpanii = "no"
params.dc_bcell = "no"
params.consolidate_epitopes = "no"
params.tcrpmhc = "no"
params.jessev = "no"
params.jessev_top_n = 3
params.include_docking_in_immunogenicity = "no"

params.allele_target_region = "Brazil"

params.netmhcpan_tcr_toprank_threshold = 0.05

// PATCH (dengue fork): defaults for previously-undeclared params, so the
// guards at lines ~1150 and ~1182 evaluate cleanly without warnings and
// the BLAST/conserved-epitope chain stays gated off by default.
params.score_conserved = "no"
params.score_t_against_b = "no"

params.b_cell_antigen_templates="P29990,P17763"  // PATCH (dengue fork): default to DENV-2 + DENV-1 polyproteins instead of alphavirus accessions

process PROCESSINPUTFASTA {
    debug true

    publishDir "${params.epitope_output_folder}/processed_fastas"

    input:
    path protein_fasta

    output:
    path("${protein_fasta.getBaseName()}_no_descriptors.fasta"), emit: fasta_wo_descriptions

    script:
    """
    #!/usr/bin/python3
    # PATCH (dengue fork): rewritten in stdlib only. Original used pandas
    # and Biopython inline which crashed with np.bool AttributeError when
    # apptainer's host-home auto-mount overlaid the container's numpy
    # with the host's numpy>=1.24. Stdlib avoids the entire issue and
    # keeps this process container-free.
    import os
    import sys

    sequences = []
    clean_protein_ids = []
    protein_descriptions = []

    cur_id = None
    cur_desc = None
    cur_seq = []
    with open("${protein_fasta}") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if cur_id is not None:
                    sequences.append(''.join(cur_seq))
                    clean_protein_ids.append(cur_id)
                    protein_descriptions.append(cur_desc)
                header = line[1:]
                parts = header.split(None, 1)
                cur_id = parts[0]
                cur_desc = header
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id is not None:
            sequences.append(''.join(cur_seq))
            clean_protein_ids.append(cur_id)
            protein_descriptions.append(cur_desc)

    out_fasta = "${protein_fasta.getBaseName()}_no_descriptors.fasta"
    with open(out_fasta, "w") as out:
        for seq, sid in zip(sequences, clean_protein_ids):
            out.write(">" + sid + "\\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\\n")

    join_dir = "${params.epitope_output_folder}/join_tables"
    if not os.path.exists(join_dir):
        os.makedirs(join_dir)

    with open(join_dir + "/input_fasta_table.tsv", "w") as out:
        out.write("\\tsequence\\tprotein_id\\tprotein_description\\n")
        for i, (seq, sid, desc) in enumerate(zip(sequences, clean_protein_ids, protein_descriptions)):
            out.write(str(i) + "\\t" + seq + "\\t" + sid + "\\t" + desc + "\\n")
    """

}

process FORMATALLELEFREQUENCIES {
    debug true

    publishDir "${params.epitope_output_folder}/join_tables"

    input:
    val regions

    output:
    path("${regions.replaceAll(/,/,"_")}_regional_allele_frequencies.tsv"), emit: regional_allele_frequencies_tsv

    script:
    """
    #!/usr/bin/python3
    # PATCH (dengue fork): rewritten in stdlib only. Original used pandas
    # which crashed in this same host-Python with the np.bool issue.
    # Replicates the same join + groupby + prevalence calculation using
    # csv module + collections.defaultdict.
    import csv
    from collections import defaultdict

    population_terms = [t.strip().lower() for t in "${regions}".split(",")]

    rows_by_allele_type = defaultdict(list)
    seen_keys = set()

    with open("${params.immunoinformatics_allele_table_path}", newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\\t')
        for row in reader:
            ds = (row.get("dataset_name") or "").lower()
            if not any(p in ds for p in population_terms):
                continue
            row_key = tuple(sorted((k, row[k]) for k in row))
            if row_key in seen_keys:
                continue
            seen_keys.add(row_key)
            allele_type = row.get("allele_type") or ""
            locus = row.get("locus") or ""
            if not allele_type or not locus or ":" not in locus:
                continue
            rows_by_allele_type[allele_type].append(locus)

    out_rows = []
    region_label = "_".join("${regions}".split(","))
    for allele_type, loci in rows_by_allele_type.items():
        # truncate locus to first two colon-separated fields
        truncated = [":".join(l.split(":")[:2]) for l in loci]
        denom = len(truncated)
        counts = defaultdict(int)
        for l in truncated:
            counts[l] += 1
        for locus, count in counts.items():
            prevalence = count / denom if denom else 0.0
            if allele_type in ("A", "B", "C"):
                locus_join = "HLA-" + locus.replace("*", "")
            else:
                locus_join = locus.replace("*", "_").replace(":", "")
            out_rows.append({
                "locus": locus,
                "count": count,
                "prevalence": prevalence,
                "locus_join": locus_join,
                "region": region_label,
            })

    out_path = "${regions.replaceAll(/,/,"_")}_regional_allele_frequencies.tsv"
    fieldnames = ["", "locus", "count", "prevalence", "locus_join", "region"]
    with open(out_path, "w", newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\\t')
        writer.writeheader()
        for i, r in enumerate(out_rows):
            r[""] = i
            writer.writerow(r)
    """

}

process CDHIT {
    debug true
    label 'CD_HIT'    

    input:
    path protein_fasta
    val similarity_threshold
    
    /*var outfile_name = "cdhit_out_sim_${similarity_threshold_str}.txt"*/

    output:
    path("cdhit_output.${similarity_threshold}"), emit: clstr_fasta
    path("cdhit_output.${similarity_threshold}.clstr"), emit: clstr_file 
    
    script:
    """
    cd-hit \
    -i $protein_fasta \
    -o cdhit_output.$similarity_threshold \
    -c $similarity_threshold \
    -l 5 \
    -T 12 \
    -g 1 \
    -aS 0.9 \
    -uS 0.1 \
    -d 200
    """   
}

process CDHITTOTSV {
    debug true
    publishDir "${params.epitope_output_folder}/cdhit_output"    

    input:    
    path clstr_file
    val similarity_threshold
    val clustering_type

    output:
    path("cdhit_out_${clustering_type}_sim_${similarity_threshold}.tsv"), emit: clstr_tsv

    script:
    """
    python3 /cdhit_clstr_to_df.py \
    $clstr_file \
    cdhit_out_${clustering_type}_sim_${similarity_threshold}.tsv
    """

}

process EPIDOPE {
    debug true 
    label 'EPIDOPE'
    publishDir "${params.epitope_output_folder}/epidope_output"

    input:
    path protein_fasta

    output:
    path("epidope_output_${protein_fasta}"), emit: epidope_output

    script:
    """
    epidope \
    -p 12 \
    -i $protein_fasta \
    -o epidope_output_$protein_fasta
    """    

}

process BEPIPRED {
    debug true

    publishDir "${params.epitope_output_folder}/bepipred_raw_output"

    input:
    path protein_fasta

    output:
    path("bepipred_output"), emit: bepipred_output

    script:
    """
    /bin/bash -c "sed '/^[^>]/s/[B|X|J|U]//g' $protein_fasta > cleaned.fasta \
    && python /bcell_standalone/predict_antibody_epitope.py -m Bepipred -f cleaned.fasta > bepipred_output"
    """
}

process BEPIPREDTOTSV {
    debug true
    publishDir "${params.epitope_output_folder}/bepipred_output"

    input:
    path bepipred_output

    output:
    path("bepipred_out.tsv"), emit: bepipred_tsv

    script:
    """
    echo Bepipred Output Path
    echo $bepipred_output

    python3 /bepipred_to_df.py \
    $bepipred_output \
    bepipred_out.tsv
    """
}

process BEPIPREDTOTSV2 {
    debug true
    publishDir "${params.epitope_output_folder}/bepipred_output"

    input:
    path bepipred_output

    output:
    path("bepipred_out.tsv"), emit: bepipred_tsv

    script:
    """
    #!/usr/bin/python3
    import sys
    import pandas as pd
    
    
    bepipred_inpath = "${bepipred_output}"
    bepipred_df_outpath = "bepipred_out.tsv"


    with open(bepipred_inpath, 'r') as f:
        rawdat = f.read().split('input: ')[1:]

        dat = [i.split('Position\\tResidue')[0] for i in rawdat]
        try:
            protname,tabledat = dat[0].split('Predicted peptides\\n')
        except:
            pass

        dfs = []

        for datchunk in dat:
            try:
                protname,tabledat = datchunk.split('Predicted peptides\\n')
                protname = protname.strip('\\n')
                dfrows = [i.split('\\t') for i in tabledat.split('\\n')]
                df = pd.DataFrame(dfrows[1:], columns=dfrows[0])
                df['protein_full'] = protname
                df['protein_id'] = protname.split(' ')[0]
                df['protein_short'] = protname.split(' ')[0].split('.')[0]
                dfs.append(df)
            except:
                pass

    all_bepipred_out = pd.concat(dfs)

    all_bepipred_out.to_csv(bepipred_df_outpath, sep='\\t', index=False)
    """
}

process MIGRATEAF2 {
    debug true

    """
    #!/usr/bin/python3
    import glob
    import os
    import shutil

    alphafold_files = glob.glob('/tmp/alphafold/*/ranked_0.pdb')
    print("Migrating AF2 Predictions:")
    print(",".join(alphafold_files))
    for alphafold_file in alphafold_files:
        protname = alphafold_file.split('/')[-2]
        outfile = f"{protname}_ranked_0.pdb"
        alphafold_dest_folder = os.path.join('/home','pathinformatics','epitope_outputs','alphafold_predictions',protname)
        alphafold_dest_file = os.path.join(alphafold_dest_folder,outfile)
        if not os.path.exists(alphafold_dest_folder):
            os.makedirs(alphafold_dest_folder)
        shutil.copy(alphafold_file, alphafold_dest_file)
    """
}

process DISCOTOPE {
    debug true

    input:
    path pdb_file

    output:
    path("${pdb_file}_discotope_out.txt"), emit: discotope_out

    script:
    """
    /discotope-1.1/discotope -f $pdb_file -chain A > ${pdb_file}_discotope_out.txt
    """
}

process DISCOTOPETOTSV {
    debug true
    publishDir "${params.epitope_output_folder}/discotope_output"

    input:
    path discotope_file

    output:
    path("${discotope_file}.tsv"), emit: discotope_out_tsv

    script:
    """
    python3 /discotope_to_tsv.py \
    $discotope_file ${discotope_file}.tsv $discotope_file
    """
}

process NETMHCPANIIEDB {
    debug true
    publishDir "${params.epitope_output_folder}/netmhcpan_i_output"

    input:
    path protein_fasta
    val allele

    output:
    path("netmhcpan_i_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_out.tsv"), emit: netmhcpan_i_tsv

    script:
    """
    /bin/bash -c "sed '/^[^>]/s/[B|X|J|U]//g' $protein_fasta > ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta \
    && /mhc_i/src/predict_binding.py netmhcpan_ba \
    $allele 9 \
    ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta > \
    netmhcpan_i_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_out.tsv"
    """
}

process NETMHCPANI {
    debug true
    publishDir "${params.epitope_output_folder}/netmhcpan_i_output"

    input:
    path protein_fasta
    val allele

    output:
    path("netmhcpan_i_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_out.xls"), emit: netmhcpan_i_xls

    script:
    """
    /bin/bash -c "sed '/^[^>]/s/[B|X|J|U]//g' $protein_fasta > ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta \
    && export NETMHCpan=/netMHCpan-4.1/Linux_x86_64/ \
    && export TMPDIR=$HOME \
    && /netMHCpan-4.1/Linux_x86_64/bin/netMHCpan \
    -BA \
    -xls \
    -a $allele \
    -l 9 \
    -t 15 \
    -v \
    -f ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta \
    -xlsfile netmhcpan_i_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_out.xls"
    """
}

process NETMHCPANIIIEDB {
    debug true
    publishDir "${params.epitope_output_folder}/netmhcpan_ii_output"

    input:
    path protein_fasta
    val allele

    output:
    path("netmhcpan_ii_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_out.tsv"), emit: netmhcpan_ii_tsv

    script:
    """
    /bin/bash -c "sed '/^[^>]/s/[B|X|J|U]//g' $protein_fasta > ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta \
    && /mhc_ii/mhc_II_binding.py netmhciipan_ba \
    $allele \
    ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta > \
    netmhcpan_ii_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_out.tsv"
    """
}

process NETMHCPANII {
    debug true
    publishDir "${params.epitope_output_folder}/netmhcpan_ii_output"

    input:
    each protein_fasta
    /*path protein_fasta*/
    val allele

    output:
    path("netmhcpan_ii_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_${protein_fasta.getName()}_out.xls"), emit: netmhcpan_ii_xls

    script:
    """
    /bin/bash -c "sed '/^[^>]/s/[B|X|J|U]//g' $protein_fasta > ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta \
    && export NETMHCIIpan=/netMHCIIpan-4.1 \
    && export TMPDIR=$HOME \
    && /netMHCIIpan-4.1/netMHCIIpan \
    -BA \
    -xls \
    -a $allele \
    -length 15 \
    -filter \
    -rankF 10 \
    -v \
    -f ${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_cleaned.fasta \
    -xlsfile netmhcpan_ii_${allele.replaceAll(/-/, "_").replaceAll(/:/,"_").replaceAll(/\*/,"_")}_${protein_fasta.getName()}_out.xls"
    """
}

process SPLITFASTAS {
    debug true
    publishDir "${params.epitope_output_folder}/split_fastas"

    input:
    path protein_fasta

    output:
    path ('split_fastas'), emit: split_fastas

    script:
    """
    python3 /split_fastas.py \
    $protein_fasta \
    /home/pathinformatics/epitope_outputs/split_fastas
    """    
}

process CONSOLIDATEEPITOPES {
    debug true

    output:
    val 1    

    script:
    """
    #!/usr/bin/python3
    import glob
    import os
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from multiprocessing import Pool

    epitope_paths = glob.glob("${params.epitope_output_folder}/*")
    # Consolidate from Epidope
    print("Consolidating Epidope")
    epidope_output_files = glob.glob("${params.epitope_output_folder}/epidope_output/epidope_output*/predicted_epitopes.csv")
    epidope_df = pd.concat([pd.read_csv(file, sep='\t') for file in epidope_output_files])[["#Gene_ID","sequence","score"]].drop_duplicates()
    ##[["#Gene_ID","sequence"]]
    epidope_df.rename(columns={"#Gene_ID":"protein_id","score":"epidope_score"}, inplace=True)
    epidope_df["type"] = "epidope"
    print(f"Epidope Table Size: {epidope_df.shape}")
    
    # Consolidate from Bepipred
    print("Consolidating Bepipred")
    bepipred_output_files = glob.glob("${params.epitope_output_folder}/bepipred_output/*.tsv")
    bepipred_df = pd.concat([pd.read_csv(file, sep='\t') for file in bepipred_output_files])[["Peptipe","Length","protein_id"]].drop_duplicates()
    ##[["Peptide","protein_id"]]
    bepipred_df.rename(columns={"Peptipe":"sequence"}, inplace=True)
    bepipred_df["type"] = "bepipred"
    bepipred_df["bepipred_score"] = 1
    # Filter Bepipred predicted epitopes by target size of 5-22 AA
    bepipred_df = bepipred_df[(bepipred_df["Length"]>=5) & (bepipred_df["Length"]<=22)]
    bepipred_df.drop("Length", axis=1, inplace=True)
    print(f"Bepipred Table Size: {bepipred_df.shape}")
    
    # Consolidate from Discotope
    print("Consolidating Discotope")
    discotope_files = glob.glob("${params.epitope_output_folder}/discotope_output/*.tsv")
    discotope_df = pd.concat([pd.read_csv(file, sep='\t') for file in discotope_files])[["mean_score","sequence","protein_name"]].drop_duplicates()
    discotope_df["protein_name"] = [i.split("_ranked")[0] for i in discotope_df["protein_name"].values]
    ##[["sequence","protein_name"]]
    discotope_df.rename(columns={"protein_name":"protein_id"}, inplace=True)
    discotope_df["type"] = "discotope"
    print(f"Discotope Table Size: {discotope_df.shape}")

    # Consolidate from NetMHCIPan (PATCH dengue fork: IEDB-wrapped NetMHCpan
    # outputs TSV not XLS; columns: allele, seq_num, start, end, length,
    # peptide, ic50, percentile_rank).
    TOPRANK = $params.netmhcpan_tcr_toprank_threshold
    toprank = str(TOPRANK).replace(".","_")
    print("Consolidating NetMHCPAN I (IEDB-format TSV)")
    netmhcpan_i_files = glob.glob("${params.epitope_output_folder}/netmhcpan_i_output/*.tsv")
    netmhcpan_i_files = [f for f in netmhcpan_i_files if 'top_' not in f and 'peptides_protein_ids' not in f]
    netmhcpan_i_dfs = []
    for file in netmhcpan_i_files:
        df = pd.read_csv(file, sep='\t')
        if df.empty:
            continue
        netmhcpan_i_dfs.append(df)
    if not netmhcpan_i_dfs:
        netmhcpan_i_df = pd.DataFrame(columns=['allele','seq_num','peptide','ic50','percentile_rank'])
    else:
        netmhcpan_i_df = pd.concat(netmhcpan_i_dfs)
    netmhcpan_i_df.rename(columns={"seq_num":"protein_id","peptide":"sequence","ic50":"netmhcpan_i_ic50","percentile_rank":"netmhcpan_i_ba_rank"}, inplace=True)
    # Synthesize columns the rest of the workflow expects
    if "netmhcpan_i_el_score" not in netmhcpan_i_df.columns:
        netmhcpan_i_df["netmhcpan_i_el_score"] = float('nan')
        netmhcpan_i_df["netmhcpan_i_el_rank"] = float('nan')
        netmhcpan_i_df["netmhcpan_i_ba_score"] = 1.0 / (netmhcpan_i_df.get("netmhcpan_i_ic50", 1.0) + 1.0)
        netmhcpan_i_df["netmhcpan_i_ave"] = netmhcpan_i_df["netmhcpan_i_ba_rank"]
        netmhcpan_i_df["netmhcpan_i_nb"] = (netmhcpan_i_df["netmhcpan_i_ba_rank"] <= 2.0).astype(int)
    netmhcpan_i_df = netmhcpan_i_df.drop_duplicates(subset=['sequence','protein_id','allele'])
    netmhcpan_i_df[["sequence", "protein_id"]].drop_duplicates().to_csv("${params.epitope_output_folder}/join_tables/netmhcpan_i_peptides_protein_ids.tsv", sep='\t')
    netmhcpan_i_df["type"] = "netmhcpan_i"
    netmhcpan_i_df.reset_index().to_feather("${params.epitope_output_folder}/netmhcpan_i_output/all_netmhcpan_i.feather")
    print(f"NetMHCPAN I Table Size: {netmhcpan_i_df.shape}")
    netmhcpan_i_toprank = netmhcpan_i_df[netmhcpan_i_df["netmhcpan_i_ba_rank"]<=TOPRANK]
    netmhcpan_i_toprank.to_csv(f"${params.epitope_output_folder}/netmhcpan_i_output/netmhcpan_i_top_{toprank}.tsv", sep='\t', index=False)
    netmhcpan_i_toprank_records = [SeqRecord(seq=Seq(s), id=s, description="") for s in list(netmhcpan_i_toprank["sequence"].unique())]
    with open(f"${params.epitope_output_folder}/netmhcpan_i_output/netmhcpan_i_top_{toprank}.fasta", "w") as netmhcpan_i_toprank_output_file:
        SeqIO.write(netmhcpan_i_toprank_records, netmhcpan_i_toprank_output_file, "fasta")


    # Consolidate from NetMHCIIPan (PATCH dengue fork: IEDB TSV format)
    TOPRANK = $params.netmhcpan_tcr_toprank_threshold
    toprank = str(TOPRANK).replace(".","_")
    print("Consolidating NetMHCPAN II (IEDB-format TSV)")
    netmhcpan_ii_files = glob.glob("${params.epitope_output_folder}/netmhcpan_ii_output/*.tsv")
    netmhcpan_ii_files = [f for f in netmhcpan_ii_files if 'top_' not in f and 'peptides_protein_ids' not in f]
    netmhcpan_ii_dfs = []
    for file in netmhcpan_ii_files:
        df = pd.read_csv(file, sep='\t')
        if df.empty:
            continue
        netmhcpan_ii_dfs.append(df)
    if not netmhcpan_ii_dfs:
        netmhcpan_ii_df = pd.DataFrame(columns=['allele','seq_num','peptide','ic50','percentile_rank'])
    else:
        netmhcpan_ii_df = pd.concat(netmhcpan_ii_dfs)
    netmhcpan_ii_df.rename(columns={"seq_num":"protein_id","peptide":"sequence","ic50":"netmhcpan_ii_ic50","percentile_rank":"netmhcpan_ii_ba_rank"}, inplace=True)
    netmhcpan_ii_df["netmhcpan_ii_el_score"] = float('nan')
    netmhcpan_ii_df["netmhcpan_ii_el_rank"] = float('nan')
    netmhcpan_ii_df["netmhcpan_ii_ba_score"] = 1.0 / (netmhcpan_ii_df.get("netmhcpan_ii_ic50", 1.0) + 1.0)
    netmhcpan_ii_df["netmhcpan_ii_ave"] = netmhcpan_ii_df["netmhcpan_ii_ba_rank"]
    netmhcpan_ii_df["netmhcpan_ii_nb"] = (netmhcpan_ii_df["netmhcpan_ii_ba_rank"] <= 5.0).astype(int)
    netmhcpan_ii_df = netmhcpan_ii_df.drop_duplicates(subset=['sequence','protein_id','allele'])
    netmhcpan_ii_df[["sequence", "protein_id"]].drop_duplicates().to_csv("${params.epitope_output_folder}/join_tables/netmhcpan_ii_peptides_protein_ids.tsv", sep='\t')
    netmhcpan_ii_df["type"] = "netmhcpan_ii"
    netmhcpan_ii_df.reset_index().to_feather("${params.epitope_output_folder}/netmhcpan_ii_output/all_netmhcpan_ii.feather")
    print(f"NetMHCPAN II Table Size: {netmhcpan_ii_df.shape}")
    netmhcpan_ii_toprank = netmhcpan_ii_df[netmhcpan_ii_df["netmhcpan_ii_ba_rank"]<=TOPRANK]
    netmhcpan_ii_toprank.to_csv(f"${params.epitope_output_folder}/netmhcpan_ii_output/netmhcpan_ii_top_{toprank}.tsv", sep='\t', index=False)
    netmhcpan_ii_toprank_records = [SeqRecord(seq=Seq(s), id=s, description="") for s in list(netmhcpan_ii_toprank["sequence"].unique())]
    with open(f"${params.epitope_output_folder}/netmhcpan_ii_output/netmhcpan_i_top_{toprank}.fasta", "w") as netmhcpan_ii_toprank_output_file:
        SeqIO.write(netmhcpan_ii_toprank_records, netmhcpan_ii_toprank_output_file, "fasta")

    # Consolidate Outputs
    all_concat_df_outdir = "${params.epitope_output_folder}/consolidated_outputs"
    if not os.path.exists(all_concat_df_outdir):
        os.makedirs(all_concat_df_outdir)

    
    
    # Create consolidated FASTA
    dfs = [epidope_df, bepipred_df, discotope_df, netmhcpan_i_df, netmhcpan_ii_df]    
    
    def create_fasta_and_txt_from_df(df):
        sequence_records = []
        df_type = df["type"].values[0]
        df_proc = df[(pd.notna(df["sequence"]))][["protein_id","sequence"]].drop_duplicates()
        txt_lines = list(set(df_proc["sequence"].unique()))

        pd.Series(df_proc["sequence"].unique()).apply(lambda x: sequence_records.append(SeqRecord(seq=Seq(str(x)), id=str(x), description=str(df_type)  )) )

        with open(f"${params.epitope_output_folder}/consolidated_outputs/consolidated_epitopes_{df_type}.fasta", 'w') as df_fasta_outfile:
            SeqIO.write(sequence_records, df_fasta_outfile, "fasta")

        with open(f"${params.epitope_output_folder}/consolidated_outputs/consolidated_epitopes_{df_type}.txt", 'w') as txt_outfile:
            txt_outfile.writelines("\\n".join(txt_lines))

    with Pool(len(dfs)) as p:
        p.map(create_fasta_and_txt_from_df, dfs)

    # Consolidate all dfs together
    for df in [epidope_df, bepipred_df, discotope_df, netmhcpan_i_df, netmhcpan_ii_df]:
        df_type = df["type"].values[0]
        df.reset_index().to_feather(f"${params.epitope_output_folder}/consolidated_outputs/consolidated_outputs_{df_type}.feather")
        df.to_csv(f"${params.epitope_output_folder}/consolidated_outputs/consolidated_outputs_{df_type}.tsv", sep='\t')
    """
}

process GATHEREPITOPEFASTAS {
    debug true

    input:
    val consolidated_output

    output:
    path("all_consolidated_epitopes.fasta"), emit: gathered_epitopes_fasta
    
    script:
    """
    cat ${params.epitope_output_folder}/consolidated_outputs/consolidated_epitopes*.fasta > \
    all_consolidated_epitopes.fasta    
    """

}

process PREPAREDATAFORJESSEV {
    debug true

    publishDir "${params.epitope_output_folder}/jessev_input"

    input:
    path allele_frequencies_table

    output:
    path("jessev_input.csv"), emit: jessev_input_csv

    script:
    """
    #!/usr/bin/python3
    import pandas as pd
    import os

    MAX_NUM_PROTEIN_IDS = 10
    THRESHOLD_VALUE = 05.00
    THRESHOLD_ATTRIBUTE = "netmhcpan_i_ba_rank"

    if not os.path.exists("${params.epitope_output_folder}/jessev_input"):
        os.makedirs("${params.epitope_output_folder}/jessev_input")

    netmhcpan_i_df = pd.read_feather("${params.epitope_output_folder}/netmhcpan_i_output/all_netmhcpan_i.feather")

    allele_freq_df = pd.read_csv("${allele_frequencies_table}", sep='\t')
    netmhcpan_i_df_full = netmhcpan_i_df.merge(allele_freq_df, left_on="allele", right_on="locus_join", how="inner")
    
    if "${params.include_docking_in_immunogenicity}"=="yes":
        pdb_episa_df = pd.read_csv("${params.epitope_output_folder}/pdb_episa_output/pdb_episa_output_table.tsv", sep='\t').rename(columns={"SEQUENCE":"sequence"})
        netmhcpan_i_df_full = netmhcpan_i_df_full.merge(pdb_episa_df, on="sequence", how="inner")
        print(netmhcpan_i_df_full.head())
        print("shape of epitope table with docking:", netmhcpan_i_df_full.shape)

    netmhcpan_i_df_full["weighted_immunogen"] = netmhcpan_i_df_full["netmhcpan_i_ba_score"] * netmhcpan_i_df_full["prevalence"]
    jess_ev_df = netmhcpan_i_df_full[netmhcpan_i_df_full[THRESHOLD_ATTRIBUTE]<THRESHOLD_VALUE][["sequence", "protein_id", "allele", "weighted_immunogen"]].drop_duplicates().groupby("sequence").agg({"protein_id":";".join, "allele":";".join, "weighted_immunogen":"sum"})
    # Because some of these epitopes appear in MANY proteins, we have to do some deduplication.
    # first we only keep the UNIQUE alleles from the aggregation which makes sense
    # second, we have to set some cap on the number of protein_ids we can carry forward to JessEV
    # and that number is captured in MAX_NUM_PROTEIN_IDS
    jess_ev_df["allele"] = jess_ev_df["allele"].apply(lambda x: ';'.join(list(set(x.split(';')))))
    jess_ev_df["protein_id"] = jess_ev_df["protein_id"].apply(lambda x: ';'.join( sorted(list(set(x.split(';'))))[:MAX_NUM_PROTEIN_IDS]))
    print(f"shape of JessEV input table: {jess_ev_df.shape}")
    jess_ev_df.reset_index().rename(columns={"weighted_immunogen":"immunogen", "sequence":"epitope", "protein_id":"proteins", "allele":"alleles"}).to_csv("jessev_input.csv", index=False)
    """
}

process RUNJESSEV {
    debug true

    input:
    path jessev_input_csv
    val iterations

    script:
    """
    #!/usr/bin/python3
    import pandas as pd
    import subprocess
    import glob
    import os
    import shutil
    import docker
    from docker.types import Mount

    client = docker.from_env()

    NUM_EPITOPES = 3
    MIN_SPACER_LEN = 4

    shutil.copy("${params.epitope_output_folder}/jessev_input/$jessev_input_csv", "${params.epitope_output_folder}/jessev_input/jessev_input_0.csv")
  
    for iter in range($iterations):
        print("running JessEV iteration:", iter)

        if iter > 0:
            input_csv_pre = f"${params.epitope_output_folder}/jessev_input/jessev_input_{iter-1}.csv"
            df_input_csv_pre = pd.read_csv(input_csv_pre)
            output_csv_pre = f"${params.epitope_output_folder}/jessev_input/jessev_output_{iter-1}.csv"
            df_output_csv_pre = pd.read_csv(output_csv_pre)
            df_input_csv_pre["jessev_used"] = [(i in df_output_csv_pre["vaccine"].values[0]) for i in df_input_csv_pre["epitope"]]
            input_csv_filtered = df_input_csv_pre[df_input_csv_pre["jessev_used"]==False].drop("jessev_used", axis=1)
            input_csv_filtered.to_csv(f"${params.epitope_output_folder}/jessev_input/jessev_input_{iter}.csv", index=False)     

        jessev_statement = f"${params.epitope_output_folder}/jessev_input/jessev_input_{iter}.csv ${params.epitope_output_folder}/jessev_input/jessev_output_{iter}.csv"
        
        client.containers.run("pmccaffrey6/jess_ev:latest",
            command=f"/opt/conda/envs/jessev/bin/python3 /JessEV/design.py -e 3 -s 4 {jessev_statement}",
            auto_remove=True,
            mounts=[Mount('/home/pathinformatics','//home/pathinformatics', type="bind")]
         )

    if not os.path.exists(f"${params.epitope_output_folder}/jessev_output/"):
        os.makedirs(f"${params.epitope_output_folder}/jessev_output/")

    pd.concat([pd.read_csv(i) for i in glob.glob(f"${params.epitope_output_folder}/jessev_input/jessev_output_*.csv")]).to_csv(
        f"${params.epitope_output_folder}/jessev_output/jessev_outputs_top_${iterations}.tsv", sep='\t', index=False)

    """
}

process TCRPMHC {
    debug true

    input:
    path input_fasta
    path template_fasta

    output:
    val "${params.epitope_output_folder}/tcrpmhc_output/${input_fasta}_TCR-pMHC.pdb", emit: tcrpmhc_output

    script:
    """
    echo $input_fasta && \
    echo $template_fasta && \

    cp $template_fasta /home/pathinformatics/epitope_outputs/tcrpmhc_output/template_$input_fasta && \
    cat $input_fasta >> /home/pathinformatics/epitope_outputs/tcrpmhc_output/template_$input_fasta && \
    cat /home/pathinformatics/epitope_outputs/tcrpmhc_output/template_$input_fasta && \

    /opt/conda/envs/TCRpMHCmodels/bin/tcrpmhc_models \
    /home/pathinformatics/epitope_outputs/tcrpmhc_output/template_$input_fasta \
    -n $input_fasta \
    -p $params.epitope_output_folder/tcrpmhc_output    
    """
}

process PDBEPISA {
    debug true

    input:
    val tcrpmhc_batch

    output:
    val tcrpmhc_batch 

    script:
    """
    mkdir "${file(tcrpmhc_batch.first()).getBaseName()}" && \
    cp ${tcrpmhc_batch.join(" ")} ${file(tcrpmhc_batch.first()).getBaseName()} && \
    /home/pathinformatics/run_mod "${file(tcrpmhc_batch.first()).getBaseName()}" $params.epitope_output_folder/pdb_episa_output
    """
}

process PDBEPISATOTABLE {
    debug true

    input:
    val pdb_files

    script:
    """
    #!/usr/bin/python3
    import pandas as pd
    import os
    import glob
    import xml
    import xml.etree.ElementTree as ET
    from Bio.PDB.Polypeptide import three_to_one

    existing_xml_files = glob.glob(f"${params.epitope_output_folder}/pdb_episa_output/*TCR-pMHC.xml")

    print("received pdb files:")
    print(f"${pdb_files}")

    print("xml files:")
    print(existing_xml_files)

    def parse_pisa_xml(xml_file):
        tree = ET.parse(xml_file)
        root = tree.getroot()
        elems = {i.tag:i for i in list(root)}

        interface_summary = elems['INTERFACESUMMARY']
        structures = list(interface_summary)
        residues = list(elems['RESIDUES'])

        structure_data = {}

        for structure, residue in zip(structures, residues):
            base_data = {i.tag:i for i in list(structure)}
            base_tags = [i.tag for i in list(structure)]

            for i in base_tags:
                if i.startswith("NUMBEROFRED"):
                    numresidues = base_data[i]
                    for i in list(numresidues):
                        base_data[i.tag] = int(i.text)
    
            sequence = ''.join([ three_to_one(list(i)[0].text.split(':')[-1].strip(' ').split(' ')[0]) for i in list(residue)])
            base_data['SEQUENCE'] = sequence
    
            base_keys = base_data.keys()
            for base_item in base_keys:
                if base_item.startswith('SOLVENTAREA'):
                    solvent_area = {i.tag:i.text for i in list(base_data[base_item])}
                elif base_item.startswith('SOLVATIONENERGY'):
                    solvation_energy = {i.tag:i.text for i in list(base_data[base_item])}
        
            for item in solvent_area.keys():
                newkey = f"SOLVENTAREA_{item}"
                base_data[newkey] = float(solvent_area[item])
        
            for item in solvation_energy.keys():
                newkey = f"SOLVATIONERGY_{item}"
                base_data[newkey] = float(solvation_energy[item])
        
            structure_data[structure.tag] = {k:v for k,v in base_data.items() if not type(v)==xml.etree.ElementTree.Element}
    
        df_out = pd.DataFrame(structure_data.values())
        return df_out

    xml_df = pd.concat([parse_pisa_xml(xml_file) for xml_file in existing_xml_files])
    print("XML DF Shape:", xml_df.shape)
    if not os.path.exists(f"${params.epitope_output_folder}/pdb_episa_output"):
        os.makedirs(f"${params.epitope_output_folder}/pdb_episa_output")
    xml_df["SEQUENCE_LENGTH"] = [len(i) for i in xml_df["SEQUENCE"]]
    xml_df[xml_df["SEQUENCE_LENGTH"]<10].to_csv(f"${params.epitope_output_folder}/pdb_episa_output/pdb_episa_output_table.tsv", sep='\t')
    """
}

process GETUNIPROTBYACCESSION {
    debug true

    publishDir "${params.epitope_output_folder}/b_cell_antigen_templates"

    input:
    val uniprot_accession

    output:
    path("${uniprot_accession}.fasta"), emit: b_cell_antigen_fasta

    script:
    """
    #!/usr/bin/python3
    import requests
    import json
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    import os

    uniprot_out = json.loads(requests.get(f"https://www.ebi.ac.uk/proteins/api/proteins/${uniprot_accession}").text)
    fastaout = [SeqRecord(seq=Seq(str(uniprot_out['sequence']['sequence'])), id=uniprot_out['accession'], description="")]
    
    if not os.path.exists("${params.epitope_output_folder}/b_cell_antigen_templates"):
        os.makedirs("${params.epitope_output_folder}/b_cell_antigen_templates")

    with open(f"${uniprot_accession}.fasta", "w") as uniprot_fasta_out:
        SeqIO.write(fastaout, uniprot_fasta_out, "fasta")
    """
}

process BLASTFROMFILES {
    debug true

    publishDir "${params.epitope_output_folder}/blast_database"

    input:
    val fasta_files
    path query_file

    output:
    path("blast_results.tsv"), emit: blast_results_tsv

    script:
    """
    echo "fasta files ${fasta_files}"
    rm -f ${params.epitope_output_folder}/blast_database/uniprot_all.fasta*
    cat ${fasta_files.join(" ")} >> ${params.epitope_output_folder}/blast_database/uniprot_all.fasta

    makeblastdb -in ${params.epitope_output_folder}/blast_database/uniprot_all.fasta -input_type fasta -dbtype prot -title "b_cell_antigen_db"

    blastp -query ${query_file} \
    -db ${params.epitope_output_folder}/blast_database/uniprot_all.fasta \
    -outfmt 6 -evalue 0.0005 > blast_results.tsv
    """

}

process FILTERPROTEINSBYBLAST {
    debug true

    publishDir "${params.epitope_output_folder}/filtered_blast_results"

    input:
    path blast_results_tsv
    path input_proteins_file

    output:
    path "*.fasta"

    script:
    """
    #!/usr/bin/python3
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    colnames = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
    blast_df = pd.read_csv(f"${blast_results_tsv}", sep='\t', header=None, names=colnames)
    blast_df = blast_df[blast_df['pident']>80.0]
    template_accessions = list(blast_df['sseqid'].unique())
    for template_accession in template_accessions:
        template_blast = blast_df[blast_df['sseqid']==template_accession]
        matched_proteins = list(template_blast['qseqid'].unique())        

        filtered_sequences = []

        for record in SeqIO.parse(f"${input_proteins_file}", "fasta"):
            if record.id in matched_proteins:
                filtered_sequences.append(record)

        print(f"number of blast hits for uniprot accession {template_accession}: {len(filtered_sequences)}")

        with open(f"blast_hits_with_{template_accession}.fasta", "w") as output_handle:
            SeqIO.write(filtered_sequences, output_handle, "fasta")
    """
}

process CLUSTALOMEGAMSA {
    debug true

    publishDir "${params.epitope_output_folder}/clustalomega_results"

    input:
    path filtered_fasta

    output:
    path "*.txt"

    script:
    """
    echo filtered_fasta $filtered_fasta && \
    clustalo -v -i $filtered_fasta --force --outfmt=clu --wrap 10000 -o ${filtered_fasta.getBaseName()}_aligned.txt
    """    
}

process CONSERVEDSEQSFROMCLUSTAL {
    debug true

    publishDir "${params.epitope_output_folder}/conserved_sequences_from_uniprot_templates"

    input:
    path clustal_output

    output:
    path("${clustal_output.getBaseName()}.fasta"), emit: clustal_conserved_fasta

    script:
    """
    #!/usr/bin/python3
    import os
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    df = pd.read_csv("${clustal_output}", sep='      ', skiprows=1, skipfooter=1, names=['seq_accession','sequence'], engine='python')

    with open("${clustal_output}", 'r') as f:
        annotations = f.readlines()[-1].strip('\\n').lstrip('                 ')
    
    indexes = list(range(len(annotations)))
    #df.drop(df.tail(1).index,inplace=True)
    
    fragments = []
    fragment = []
    
    for a,i in zip(annotations, indexes):
        if (a != " "):
            fragment.append(i)
        else:
            if len(fragment) >= 9:
                fragments.append(fragment)
            fragment = []
        
    sequence = df['sequence'].values[0]

    output_seqrecords = []
    
    for idx,fragment in enumerate(fragments):
        seq = sequence[fragment[0]:fragment[-1]]
        clustal_filename = "${clustal_output}"
        output_seqrecords.append(SeqRecord(seq=Seq(seq), name="", id=f"{clustal_filename.split('_')[3]}_fragment_{idx}", description=""))
        #output_seqrecords.append(SeqRecord(seq=Seq(seq), name="", id=f"fragment_{idx}", description=""))

    if not os.path.exists("${params.epitope_output_folder}/conserved_sequences_from_uniprot_templates"):
        os.makedirs("${params.epitope_output_folder}/conserved_sequences_from_uniprot_templates")

    with open("${clustal_output.getBaseName()}.fasta", "w") as output_handle:
        SeqIO.write(output_seqrecords, output_handle, "fasta")
    """
}

process COMBINEDCLUSTALEPITOPEFASTAS {
    debug true

    publishDir "${params.epitope_output_folder}/conserved_sequences_from_uniprot_templates"

    input:
    val conserved_fastas

    output:
    path("clustal_combined_output.fasta"), emit: clustal_combined_output_fasta

    script:
    """
    echo ${conserved_fastas}
    cat ${conserved_fastas.join(" ")} >> clustal_combined_output.fasta
    """
}

process COMPILEFASTAFROMBCELLWORK {
    debug true

    publishDir "${params.epitope_output_folder}/b_cell_combined_fasta"

    input:
    val b_cell_tsvs

    output:
    path("combined_b_cell_epitopes.fasta")

    script:
    """
    #!/usr/bin/python3
    import pandas as pd
    import os
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO

    input_files = "${b_cell_tsvs}"

    fasta_seqrecords = []    

    for input_file in [i.strip('[').strip(']').strip(' ') for i in input_files.split(",")]:
        print("input file:", input_file)

        if ("epidope" in input_file):
            epidope_path = os.path.join(input_file,"predicted_epitopes.csv")
            
            epidope_df = pd.read_csv(epidope_path, sep='\\t')
            for idx,row in epidope_df.iterrows():
                record = SeqRecord(seq=Seq(row['sequence']), id=row['#Gene_ID'], name="", description="")
                if (len(record.seq) >= 5):
                    fasta_seqrecords.append(record)

        elif ('bepipred' in input_file):
            print('BEPIPRED')
            bepipred_df = pd.read_csv(input_file, sep='\\t')
            #print(bepipred_df)
            for idx,row in bepipred_df[pd.notna(bepipred_df['Peptipe'])].iterrows():
                record = SeqRecord(seq=Seq(row['Peptipe']), id=row['protein_id'], name="", description="")
                if (len(record.seq) >= 5):
                    fasta_seqrecords.append(record)
        
        #print(fasta_seqrecords)

    with open("combined_b_cell_epitopes.fasta", "w") as output_handle:
        SeqIO.write(fasta_seqrecords, output_handle, "fasta")
    
    """
}

workflow {
    protein_fasta_ch = Channel.fromPath(params.protein_file)
    protein_fasta_value_ch = file(params.protein_file)
    protein_fasta_clean_ch = PROCESSINPUTFASTA(protein_fasta_value_ch)
    protein_fasta_clean_ch.view()

    tcrpmhc_templates_value_ch = file(params.tcrpmhc_template_file)

    /*COLLECT B-CELL ANTIGEN TEMPLATES*/
    // PATCH (dengue fork): the entire B-cell antigen template chain
    // (UniProt fetch + BLAST + Clustal-Omega + conserved epitopes) is
    // ONLY needed when BOTH a B-cell scorer is enabled AND
    // score_conserved == "yes" (those are the only branches that
    // consume conserved_epitopes_ch — see lines ~1182, ~1190 below).
    // Without this guard the chain fires unconditionally and requires
    // blast_latest.sif + clustalomega_latest.sif that are not built
    // by default.
    if (params.score_conserved == "yes" &&
        (params.bepipred == "yes" || params.epidope == "yes" || params.dc_bcell == "yes")) {
        b_cell_antigen_fastas = GETUNIPROTBYACCESSION(Channel.from(params.b_cell_antigen_templates.split(",")))
        blast_results_tsv_ch = BLASTFROMFILES(b_cell_antigen_fastas.collect(), protein_fasta_value_ch)
        filtered_fasta_ch = FILTERPROTEINSBYBLAST(blast_results_tsv_ch, protein_fasta_value_ch).flatten()
        clustal_omega_output_ch = CLUSTALOMEGAMSA(filtered_fasta_ch)
        clustal_formatted_output_ch = CONSERVEDSEQSFROMCLUSTAL(clustal_omega_output_ch)
        conserved_epitopes_ch = COMBINEDCLUSTALEPITOPEFASTAS(clustal_formatted_output_ch.collect())
        conserved_epitopes_ch.view()
    }

    /* CALCULATE ALLELE FREQUENCY TABLES */
    allele_frequencies_table_ch = FORMATALLELEFREQUENCIES(params.allele_target_region)
    allele_frequencies_table_ch.view()

    /*protein_fasta_splits_value_ch = Channel.fromPath(protein_fasta_value_ch).splitFasta(by: params.netmhcpan_chunk_size, file: true).collect()*/
    protein_fasta_splits_value_ch = protein_fasta_clean_ch.splitFasta(by: params.netmhcpan_chunk_size, file: true).collect()

    alphafold_pdb_ch = Channel.fromPath(params.alphafold_pdb_folder)

    mhc_i_alleles_ch = Channel.from('HLA-A01:01','HLA-A02:01','HLA-A03:01','HLA-A24:02','HLA-A26:01','HLA-B07:02','HLA-B08:01','HLA-B15:01','HLA-B27:05','HLA-B39:01','HLA-B40:01','HLA-B58:01')
    mhc_ii_alleles_ch = Channel.from('DRB1_0301','DRB1_0701','DRB1_1501','DRB3_0101','DRB3_0202','DRB4_0101','DRB5_0101')    

    /*split_fastas_ch = SPLITFASTAS(protein_fasta_ch)*/
    split_fastas_ch = Channel.fromPath('/home/pathinformatics/epitope_outputs/split_fastas/*')
    
    /* CDHIT INPUT PROTEINS */
    if (params.cdhit_input_proteins == "yes") {
        cdhit_out_ch = CDHIT(protein_fasta_ch, params.cdhit_similarity_threshold)
        CDHITTOTSV(cdhit_out_ch.clstr_file, params.cdhit_similarity_threshold, "input_proteins")
    }
    /* B-CELL SCORING */
    if (params.bepipred == "yes") {
        if (params.score_conserved == "yes") {
            bepipred_out_ch = BEPIPRED(conserved_epitopes_ch)
        } else {
            bepipred_out_ch = BEPIPRED(protein_fasta_clean_ch)
        }
        bepipred_tsv_ch = BEPIPREDTOTSV2(bepipred_out_ch.bepipred_output)
    }
    if (params.epidope == "yes") {
        if (params.score_conserved == "yes") {
            epidope_tsv_ch = EPIDOPE(conserved_epitopes_ch)
        } else {
            epidope_tsv_ch = EPIDOPE(protein_fasta_clean_ch.splitFasta(by: 500, file: true))
        }
    }
    if (params.dc_bcell == "yes") {
        MIGRATEAF2()
        discotope_out_ch = DISCOTOPE(alphafold_pdb_ch)
        DISCOTOPETOTSV(discotope_out_ch)
    }
    // PATCH (dengue fork): only invoke COMPILEFASTAFROMBCELLWORK when at
    // least one B-cell scorer ran. Upstream unconditionally referenced
    // bepipred_tsv_ch and epidope_tsv_ch which crashes when both stages
    // are disabled (e.g. partial smoke runs).
    if (params.bepipred == "yes" && params.epidope == "yes") {
        COMPILEFASTAFROMBCELLWORK(bepipred_tsv_ch.mix(epidope_tsv_ch).collect())
    } else if (params.bepipred == "yes") {
        COMPILEFASTAFROMBCELLWORK(bepipred_tsv_ch.collect())
    } else if (params.epidope == "yes") {
        COMPILEFASTAFROMBCELLWORK(epidope_tsv_ch.collect())
    }

    /* T-CELL SCORING */
    // PATCH (dengue fork): Route to IEDB-wrapped NetMHCpan (NETMHCPANIIEDB /
    // NETMHCPANIIIEDB) instead of the direct DTU NetMHCpan binaries. The DTU
    // binary tarballs are academic-license per-user and are not on the cluster.
    // The IEDB MHC tools redistribute NetMHCpan with their own license; same
    // underlying predictor, slightly different output format (TSV not XLS).
    // CONSOLIDATEEPITOPES is patched below to consume the TSV format.
    if (params.netmhcpani == "yes") {
        NETMHCPANIIEDB(protein_fasta_clean_ch, mhc_i_alleles_ch)
    }
    if (params.netmhcpanii == "yes") {
        NETMHCPANIIIEDB(protein_fasta_splits_value_ch.flatten(), mhc_ii_alleles_ch)
    }
    
    if (params.consolidate_epitopes == "yes") {
        consolidated_epitopes_output_ch = CONSOLIDATEEPITOPES()
        consolidated_epitopes_fasta_ch = GATHEREPITOPEFASTAS(consolidated_epitopes_output_ch)
        
        cdhit_epitopes_out_ch = CDHIT(consolidated_epitopes_fasta_ch, params.cdhit_similarity_threshold)
        CDHITTOTSV(cdhit_epitopes_out_ch.clstr_file, params.cdhit_similarity_threshold, "epitopes")
    }
    if (params.score_t_against_b == "yes") {
        println("SCORE T AGAINST B")
    }
    if (params.tcrpmhc == "yes") {
        protein_records_ch = Channel.fromPath("${params.epitope_output_folder}/netmhcpan_i_output/netmhcpan_i_top*.fasta").splitFasta(by: 1, file:true).take(100)
        tcrpmhc_pdbs = TCRPMHC(protein_records_ch, tcrpmhc_templates_value_ch)
        tcrpmhc_batches = tcrpmhc_pdbs.collate(3)
        tcrpmhc_batches.view()
        pdb_episa_output_ch = PDBEPISA(tcrpmhc_batches)
        PDBEPISATOTABLE(pdb_episa_output_ch.collect())
    }
    if (params.jessev == "yes") {
        jessev_input_file_ch = PREPAREDATAFORJESSEV(allele_frequencies_table_ch)
        jessev_input_file_ch.view()
        RUNJESSEV(jessev_input_file_ch, params.jessev_top_n)
    }
}
