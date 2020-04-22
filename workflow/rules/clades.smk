rule proteom_preclustering:
    input : "{base}/{taxononomy}/{taxon}.gids"
    output : representatives = "{base}/{taxononomy}/{taxon}.preclustering.faa",
             preclusters = "{base}/{taxononomy}/{taxon}.preclustering",
             full_proteom = "{base}/{taxononomy}/{taxon}.faa.gz"
    params : script = "workflow/scripts/proteom_preclustering.py"
    conda : "../envs/proteom_preclustering.yaml"
    shell : """
        python {params.script} {input} {output.full_proteom} {output.representatives} {output.preclusters} {threads}
        """

rule proteom_annotations:
    input : aas = "{base}/{taxononomy}/{taxon}.preclustering.faa",
            preclusters = "{base}/{taxononomy}/{taxon}.preclustering"
    output : emap_out = "{base}/{taxononomy}/{taxon}.emapper"
    params : script = "workflow/scripts/proteom_annotations.py"
    conda : "../envs/proteom_annotation.yaml"
    shell : """
        python {params.script} {input.aas} {input.preclusters} {output.emap_out} {threads}
        """
