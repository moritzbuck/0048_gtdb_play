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
