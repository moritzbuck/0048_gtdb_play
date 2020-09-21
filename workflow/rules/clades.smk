rule proteom_preclustering:
    input : "{base}/{taxononomy}/{taxon}.gids"
    output : representatives = "{base}/{taxononomy}/{taxon}.preclustering.faa",
             preclusters = "{base}/{taxononomy}/{taxon}.preclustering",
             full_proteom = "{base}/{taxononomy}/{taxon}.faa.gz"
    params : script = "workflow/scripts/proteom_preclustering.py"
    conda : "../envs/proteom_preclustering.yaml"
    threads : 24
    shell : """
        python {params.script} {input} {output.full_proteom} {output.representatives} {output.preclusters} {threads}
        """

rule silix_clade:
    input : gids = "{base}/{taxononomy}/{taxon}.gids",
            aas = "{base}/{taxononomy}/{taxon}.preclustering.faa",
            preclusters = "{base}/{taxononomy}/{taxon}.preclustering"
    output : silix_clusts = "{base}/{taxononomy}/{taxon}.silix.clusters"
    params : script = "workflow/scripts/silix.py"
    conda : "../envs/silix.yaml"
    threads : 24
    shell : """
        python3 {params.script} {input.gids} {input.aas} {input.preclusters} {output.silix_clusts} {threads}
        """

rule proteom_annotation:
    input : gids = "{base}/{taxononomy}/{taxon}.gids",
            aas = "{base}/{taxononomy}/{taxon}.preclustering.faa",
            preclusters = "{base}/{taxononomy}/{taxon}.preclustering"
    output : emap_out = "{base}/{taxononomy}/{taxon}.emapper"
    params : script = "workflow/scripts/proteom_annotation.py"
    conda : "../envs/proteom_annotation.yaml"
    threads : 24
    shell : """
        python3 {params.script} {input.gids} {input.aas} {input.preclusters} {output.emap_out} {threads}
        """

rule motupan_COGs:
    input :  gids = "{base}/{taxononomy}/g__{genus}/s__{species}/s__{species}.gids",
             silix_clusts = "{base}/{taxononomy}/g__{genus}/g__{genus}.silix.clusters"
    output : motupan_out = "{base}/{taxononomy}/g__{genus}/s__{species}/s__{species}.ab_initio.motupan.json"
    params : script = "workflow/scripts/motupan.py", checkm_file = "data/checkm_file.txt"
    conda : "../envs/motulizer.yaml"
    threads : 1
    shell : """
        python3 '{params.script}' '{input.gids}' '{input.silix_clusts}' '{output.motupan_out}' '{params.checkm_file}'
        """

rule motupan_emapper:
    input :  gids = "{base}/{taxononomy}/g__{genus}/s__{species}/s__{species}.gids",
             emapper_results = "{base}/{taxononomy}/g__{genus}/g__{genus}.emapper"
    output : motupan_out = "{base}/{taxononomy}/g__{genus}/s__{species}/s__{species}.emapper.motupan.json"
    params : script = "workflow/scripts/emap_motupan.py", checkm_file = "data/checkm_file.txt"
    conda : "../envs/motulizer.yaml"
    threads : 1
    shell : """
        python3 '{params.script}' '{input.gids}' '{input.emapper_results}' '{output.motupan_out}' '{params.checkm_file}'
        """
