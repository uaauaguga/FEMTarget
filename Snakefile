import json
import os

#load all genome ids
genome_ids = open(config["genome-ids"]).read().strip().split("\n")
# load query genome ids
query_ids = open(config["query-ids"]).read().strip().split("\n")
print(query_ids)
for query_id in query_ids:
    if query_id not in genome_ids:
        print(f"{query_id} not in genome set, insert it.")
        genome_ids.append(query_id)


# define input and output directory
indir = config["indir"]
outdir = config["outdir"]

sRNA_ids_considered0 = None
if "sRNA-ids" in config:
    sRNA_ids_considered0 = open(config["sRNA-ids"]).read().strip().split("\n")

genome_ids_by_srna = {}
sRNA_ids_considered = []
for genome_id in genome_ids:
    for fasta in os.listdir(f"{indir}/sRNAs/{genome_id}"):
        #print(fasta)
        sRNA_id = fasta[:fasta.rfind(".")]
        if (sRNA_ids_considered0 is not None) and (sRNA_id not in sRNA_ids_considered0):
            continue
        if sRNA_id not in genome_ids_by_srna:
            genome_ids_by_srna[sRNA_id] = []
        genome_ids_by_srna[sRNA_id].append(genome_id) 

for sRNA_id in genome_ids_by_srna:
    if len(genome_ids_by_srna[sRNA_id]) >= 4:
        #if (sRNA_ids_considered0 is not None) and (sRNA_id not in sRNA_ids_considered0):
        #    continue
        sRNA_ids_considered.append(sRNA_id)   

#print(sRNA_ids_considered)
weights = config["weights"]

params_string = []
weights_string = []

for k,v in weights.items():
    k = k.replace(" ",".")
    params_string.append(f"{k}-{v}")
    weights_string.append(f"{k}-{v}")
weights_string = "_".join(weights_string)
denoise = {}
combine = {}
denoise = config["denoise"]
for k,v in denoise.items():    
    params_string.append(f"{k}-{v}")

params_string = "_".join(params_string)

params_string = config.get("hfq","nohfq") + "--" + params_string

genome_set = config["genome-set-name"]
def get_output(wildcards):
    paths = []
    for genome_id in genome_ids:
        sRNA_ids = []
        for fa in os.listdir(f"{indir}/sRNAs/{genome_id}/"):
            if fa[:-3] not in sRNA_ids_considered:
                continue
            sRNA_ids.append(fa[:-3])
    paths += expand(outdir + "/{query_id}/homologs/{genome_set}/best.hit.txt",query_id=query_ids, genome_set = genome_set)
    paths += [outdir + f"/proteins/{genome_set}/table.txt"]     
    if "denoise" in config:
        paths += expand(outdir + f"/comparative-analysis-denoised-ml/{genome_set}/" + params_string + "/{sRNA_id}.txt", sRNA_id = sRNA_ids_considered)
    else:
        paths += expand(outdir + f"/comparative-analysis-combined/{genome_set}/" + params_string + "/{sRNA_id}.txt", sRNA_id = sRNA_ids_considered)
    paths += expand(outdir + f"/checkpoints/{genome_set}/" + params_string + "/{sRNA_id}.txt",sRNA_id = sRNA_ids_considered)
    return paths
                
rule all:
    input:
        score = get_output

rule prepare_bed:
    input:
        gff =  indir + "/gff/{genome_id}.gff"
    output:
        bed = indir + "/CDS/{genome_id}.bed"
    shell:
        """
        scripts/gff2bed.py --gff {input.gff} --bed {output.bed} --feature CDS --name Name,locus_tag,gene 
        """

rule prepare_protein:
    input:
        fasta = indir + "/fasta/{genome_id}.fa",
        bed = indir + "/CDS/{genome_id}.bed"
    output:
        protein = outdir + "/{genome_id}/proteins.faa",
    shell:
        """
        scripts/extract-protein-sequence-from-genome.py -i {input.bed} -g {input.fasta} -o {output.protein}
        """        

rule extract_leader:
    input:
        fasta = indir +"/fasta/{genome_id}.fa",
        CDS = indir + "/CDS/{genome_id}.bed"
    output:
        leader = outdir + "/{genome_id}/leader.fa",
        fai = indir +"/fasta/{genome_id}.fa.fai"
    shell:
        """
        scripts/extract-leader-sequences.py -i {input.CDS} -g {input.fasta} -o {output.leader} 
        """ 

rule energy_scoring:
    input:
        sRNA = indir + "/sRNAs/{genome_id}/{sRNA_id}.fa",
        leaders = outdir + "/{genome_id}/leader.fa"
    output:
        energy = outdir + "/{genome_id}/energies/{sRNA_id}.txt"
    log:
        log = outdir + "/{genome_id}/energies/{sRNA_id}.log"
    threads: 1 if config["energy-scorer"] == "RNAplex" else 10
    params:
        method = config["energy-scorer"]
    shell:
        """
        scripts/calculate-energy-score.py -rs "{input.sRNA}" -ts "{input.leaders}" -o "{output.energy}" -j {threads} -m {params.method} > {log.log} 2>&1
        """


rule energy_normalization:
    input:
        energy = outdir + "/{genome_id}/energies/{sRNA_id}.txt",
    output:
        energy = outdir + "/{genome_id}/energies.normalized/{sRNA_id}.txt"
    shell:
        """
        scripts/normalize-energy.py -i "{input.energy}" -o "{output.energy}" 
        """


rule combine_scores_wthfq:
    input:
        energy = outdir + "/{genome_id}/energies.normalized/{sRNA_id}.txt",
        hfq = indir + "/" + config["hfq"] + "/{genome_id}.txt"
    output:
        score = outdir + "/{genome_id}/combined/{sRNA_id}.txt"
    log:
        log = outdir + "/{genome_id}/combined/{sRNA_id}.log"
    shell:
        """   
        scripts/combine-scores.py -es "{input.energy}" -hs {input.hfq} -o "{output.score}" > {log.log} 2>&1
        """

rule combine_scores_wohfq:
    input:
        energy = outdir + "/{genome_id}/energies.normalized/{sRNA_id}.txt",
    output:
        score = outdir + "/{genome_id}/combined.wo.hfq/{sRNA_id}.txt"
    log:
        log = outdir + "/{genome_id}/combined.wo.hfq/{sRNA_id}.log"
    shell:
        """   
        scripts/combine-scores.py -es "{input.energy}" -o "{output.score}" > "{log.log}" 2>&1
        """

rule combine_proteins:
    input:
        proteins = expand(outdir + "/{genome_id}/proteins.faa",genome_id=genome_ids)
    output:
        proteins = outdir + "/proteins/{genome_set}/db.faa"
    params:
        od = outdir + "/proteins"
    run:
        import os
        print("combine proteins ...")
        fout = open(output.proteins,"w")
        for path in input.proteins:
            genome_id = path.split("/")[-2]
            with open(path) as f:
                for line in f:
                    if line.startswith(">"):
                        line = ">" + genome_id + ":" + line[1:]            
                    fout.write(line)
        fout.close()

rule build_mmseqs_db:
    input:
        proteins = outdir + "/proteins/{genome_set}/db.faa"
    output:
        db = outdir + "/proteins/{genome_set}/db"
    shell:
        """
        export PATH=$PATH:/BioII/lulab_b/jinyunfan/miniforge3/envs/rna-analysis/bin
        mmseqs createdb --dbtype 1 {input.proteins} {output.db}
        """

rule homolog_search:
    input:
        proteins = outdir + "/proteins/{genome_set}/db.faa",
        db = outdir + "/proteins/{genome_set}/db",
        query = outdir + "/{query_id}/proteins.faa" 
    output:
        hits = outdir + "/{query_id}/homologs/{genome_set}/best.hit.txt",
    log:
        log = outdir + "/{query_id}/homologs/{genome_set}.log"
    params:
        od = outdir + "/{query_id}/homologs/{genome_set}",
        db = outdir + "/proteins/{genome_set}/db"
    threads: 5
    shell:
        """
        export PATH=$PATH:/BioII/lulab_b/jinyunfan/miniforge3/envs/rna-analysis/bin
        scripts/protein-homolog-search.py -q {input.query} -db {params.db} -p {input.proteins} -od {params.od} > {log.log} 2>&1
        """


rule collapse:
    input:
        hits = expand(outdir + "/{query_id}/homologs/" + f"{genome_set}/best.hit.txt",query_id=query_ids),
        query_ids = config["query-ids"]
    output:
        table = outdir + f"/proteins/{genome_set}/table.txt",
        collapse = outdir + f"/proteins/{genome_set}/collapse.txt"
    params:
        wd = outdir,
        genome_set = genome_set
    shell:
        """
        scripts/merge-protein-homolog-hits.py -wd {params.wd} -qi {input.query_ids} -o {output.table} \
        -c {output.collapse} --genome-set {params.genome_set}
        """

withhfq = False
if ("hfq leader zscore" in weights) and (weights["hfq leader zscore"] > 0):
    withhfq = True
       
def get_score_requirements(wildcards):
    paths = []
    sRNA_id = wildcards.sRNA_id
    for g in genome_ids_by_srna[wildcards.sRNA_id]:
        if withhfq:
            paths.append(outdir + f"/{g}/combined/{sRNA_id}.txt")
        else:
            paths.append(outdir + f"/{g}/combined.wo.hfq/{sRNA_id}.txt")
    return paths

rule collate_scores_by_homolog:
    input:
        srna = indir + "/hits.by.sRNAs/{sRNA_id}.fa",
        score = get_score_requirements, 
        table = outdir + f"/proteins/{genome_set}/table.txt",
        genome_ids = config["genome-ids"]
    output:
        scores = outdir + "/comparative-analysis/{genome_set}/{sRNA_id}.txt" if withhfq else outdir + "/comparative-analysis.wo.hfq/{genome_set}/{sRNA_id}.txt"
    params:
        wd = outdir,
        genome_set = genome_set,
        no_hfq = ["","--no-hfq"][int(not withhfq)]
    shell:
        """
        scripts/collate-homologous-pairs.py -pg {input.table} --srna "{input.srna}" \
        --input-directory {params.wd} --output "{output.scores}" {params.no_hfq} -gi {input.genome_ids}
        """ 

rule prepare_weights:
    input:
    output:
        weights = outdir + "/weights." + weights_string + ".json"
    params:
        weights = json.dumps(weights)
    run:
        with open(output.weights,"w") as f:
            f.write(params.weights)

rule extract_16S_rRNA_sequences:
    input:
        fasta = indir +"/fasta/{genome_id}.fa"
    output:
        fasta = indir + "/SSU-rRNA/" + "{genome_id}.fa"
    log:
        log = indir + "/SSU-rRNA/" + "{genome_id}.log"
    shell:
        """
        export PATH=$PATH:/BioII/lulab_b/jinyunfan/miniforge3/envs/rna-analysis/bin
        scripts/extract-16S-rRNAs.py --fasta {input.fasta} --output {output.fasta} > {log.log} 2>&1
        """

rule build_16S_rRNA_tree:
    input:
        fastas = expand(indir + "/SSU-rRNA/{genome_id}.fa", genome_id=genome_ids),
        genome_ids = config["genome-ids"]
    output:
        tree = outdir + f"/phylogeny/{genome_set}/16S-rRNA.nwk",
        msa = outdir + f"/phylogeny/{genome_set}/16S-rRNA.afa"
    params:
        indir = indir + "/SSU-rRNA",
        outdir = outdir + f"/phylogeny/{genome_set}"
    log:
        log =  outdir + f"/phylogeny/{genome_set}.log"
    shell:
        """
        export PATH=$PATH:/BioII/lulab_b/jinyunfan/miniforge3/envs/rna-analysis/bin
        scripts/build-16S-rRNA-tree.py -id {params.indir} -od {params.outdir} -gi {input.genome_ids} > {log.log} 2>&1
        """           

rule extract_representative_genomes:
    input:
        msa = outdir + f"/phylogeny/{genome_set}/16S-rRNA.afa"
    output:
        all_fasta = outdir + f"/phylogeny/{genome_set}/16S-rRNA.fa",
        rep_fasta = outdir + f"/phylogeny/{genome_set}/16S-rRNA.rep.fa",
        rep = outdir + f"/phylogeny/{genome_set}/representative-genome-ids.txt"
    shell:
        """
        export PATH=$PATH:/BioII/lulab_b/jinyunfan/miniforge3/envs/rna-analysis/bin
        esl-reformat --informat afa fasta {input.msa} > {output.all_fasta}
        cd-hit-est -i {output.all_fasta} -o {output.rep_fasta} -c 0.9999999 -r 0 -d 1000
        cat {output.rep_fasta} | grep '>' | sed 's/>//g' > {output.rep}
        """

rule denosing:
    input:
        scores =  outdir + "/comparative-analysis/{genome_set}/{sRNA_id}.txt" if withhfq else outdir + "/comparative-analysis.wo.hfq/{genome_set}/{sRNA_id}.txt",
        tree = outdir + f"/phylogeny/{genome_set}/16S-rRNA.nwk",
        weights = outdir + "/weights." + weights_string + ".json"
    output:
        scores = outdir + "/comparative-analysis-denoised-ml/{genome_set}/" + params_string + "/{sRNA_id}.txt"
    log:
        log = outdir + "/comparative-analysis-denoised-ml/{genome_set}/" + params_string + "/{sRNA_id}.log"
    threads: 10
    params:
        srm = denoise.get("srm",0.08),
        nvm = denoise.get("nvm",0.004),
        normalize = {0:"",1:"--normalize"}[int(config.get("normalize",0))],
        rescale = ["","--scale-tree"][denoise.get("rescale",0)],
        reroot = ["","--reroot"][denoise.get("reroot",0)]
    shell:
        """
        export PATH=/BioII/lulab_b/jinyunfan/miniforge3/envs/rna-analysis/bin:$PATH
        scripts/comparative-scoring-ml.py --input {input.scores} --output {output.scores} --tree  {input.tree} -sr {params.srm} -nv {params.nvm}  --jobs 10 -w {input.weights} --reroot  > {log.log} 2>&1
        """

rule combine_interaction_scores:
    input:
        scores =  outdir + "/comparative-analysis/{genome_set}/{sRNA_id}.txt" if withhfq else outdir + "/comparative-analysis.wo.hfq/{genome_set}/{sRNA_id}.txt",
        tree = outdir + f"/phylogeny/{genome_set}/16S-rRNA.nwk",
        weights = outdir + "/weights." + weights_string + ".json"
    output:
        scores = outdir + "/comparative-analysis-combined/{genome_set}/" + params_string + "/{sRNA_id}.txt"
    params:
        correlation = combine.get("correlation",0)
    shell:
        """
        scripts/aggregate-scores.py -i "{input.scores}" -o "{output.scores}" -t {input.tree} -w {input.weights} -c {params.correlation}
        """

rule extract_score:
    input:
        comparative_scores = outdir + f"/comparative-analysis-denoised-ml/{genome_set}/" + params_string + "/{sRNA_id}.txt" if "denoise" in config else outdir + f"/comparative-analysis-combined/{genome_set}/" + params_string + "/{sRNA_id}.txt",
        genome_scores = outdir + f"/comparative-analysis/{genome_set}" + "/{sRNA_id}.txt" if withhfq else outdir + f"/comparative-analysis.wo.hfq/{genome_set}" + "/{sRNA_id}.txt",        
        collapse =  outdir + f"/proteins/{genome_set}/collapse.txt",
        weights = outdir + "/weights." + weights_string + ".json",
        rep = outdir + f"/phylogeny/{genome_set}/representative-genome-ids.txt" 
    output:
        checkpoint = outdir + f"/checkpoints/{genome_set}/" + params_string + "/{sRNA_id}.txt"
    params:
        wd = outdir,
        name = "{sRNA_id}",
        tag = genome_set + "." + params_string,
        conservation_weight = config.get("conservation",0),
        normalize = {0:"",1:"--normalize"}[int(config.get("normalize",0))]
    log:
        log = outdir + f"/checkpoints/{genome_set}/" + params_string + "/{sRNA_id}.log"
    shell:
        """
        scripts/extract-final-score.py --comparative-scores "{input.comparative_scores}" --genome-scores "{input.genome_scores}" -cw {params.conservation_weight} \
        --collapse {input.collapse} -w {input.weights} --working-directory {params.wd}  {params.normalize} \
        --srna-name "{params.name}" --tag {params.tag}".ml" -rgi {input.rep} > "{log.log}" 2>&1 && touch "{output.checkpoint}"
        """
