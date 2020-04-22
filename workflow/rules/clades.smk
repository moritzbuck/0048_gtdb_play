rule proteom_preclustering:
    input : "{base}/{taxononomy}/{taxon}.gids":
    output : representatives = "{base}/{taxononomy}/{taxon}.preclustering.faa",
             preclusters = "{base}/{taxononomy}/{taxon}.preclustering",
             full_proteom = "{base}/{taxononomy}/{taxon}.faa.gz",
    params : script = "workflow/scripts/proteom_preclustering.py",
    env : "workflow/envs/proteom_preclustering.yaml"
    shell : """
        python {script} {input} {output.full_proteom} {output.representatives} {output.preclusters} {threads}
        """
