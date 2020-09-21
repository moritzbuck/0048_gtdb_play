rule analysis_bmft:
    input : db = "{base}/dbs/{tax_level}/{traits}/database.json"
    output : bmft = "{base}/analyses/{tax_level}/{traits}/master_table.csv"
    params : script = "workflow/scripts/bmft.py"
    conda : "../envs/python_std.yaml"
    threads : 1
    shell : """
        python3 {params.script} {input.db} {output.bmft} {wildcards.tax_level} {wildcards.traits}
        """


rule anot_matrices:
    input : db = "{base}/dbs/{tax_level}/{traits}/database.json",
            annot_db = "{base}/dbs/{tax_level}/annotation_databases/cogcat-{traits}.json"
    output : ec_table = "{base}/analyses/{tax_level}/{traits}/annotations/ec.csv",
             cogcat_table = "{base}/analyses/{tax_level}/{traits}/annotations/cogcat.csv",
             cazy_table = "{base}/analyses/{tax_level}/{traits}/annotations/cazy.csv",
             ko_table = "{base}/analyses/{tax_level}/{traits}/annotations/ko.csv",
    params : script = "workflow/scripts/annot_matrices.py"
    conda : "../envs/python_std.yaml"
    threads : 1
    shell : """
        python3 {params.script} {input.db} {input.annot_db} {output.cogcat_table} {wildcards.tax_level} {wildcards.traits}
        """




rule figs:
    input : bmft = "{base}/analyses/{tax_level}/{traits}/master_table.csv",
            pairs = "{base}/analyses/{tax_level}/{traits}/pairs_table.csv",
            cogcats = "{base}/analyses/{tax_level}/{traits}/annotations/cogcat.csv"
    output : fig1 = "{base}/analyses/{tax_level}/{traits}/plots/Fig1.pdf"
    params : script = "workflow/scripts/figs.R"
    conda : "../envs/R_std.yaml"
    threads : 1
    shell : """
        Rscript {params.script} {input.bmft} {input.pairs} {input.cogcats} {output.fig1} {wildcards.tax_level} {wildcards.traits}
        """

rule gtdbtk_tree:
    input : db = "{base}/dbs/{tax_level}/{traits}/database.json"
    output : arc_tree = "{base}/analyses/{tax_level}/{traits}/tree/gtdbtk_denovo_arc_rep.nwk",
             bac_tree = "{base}/analyses/{tax_level}/{traits}/tree/gtdbtk_denovo_bac_rep.nwk",
             dist_json = "{base}/analyses/{tax_level}/{traits}/tree/rep_dist.json"
    params : script = "workflow/scripts/make_gtdtk_tree.py"
    conda : "../envs/gtdbtk.yaml"
    threads : 24
    shell : """
        python3 {params.script} {input.db} {output.arc_tree} {output.bac_tree} {output.dist_json} {wildcards.tax_level} {wildcards.traits}
        """

rule analysis_pairwise:
    input : bmft = "{base}/analyses/{tax_level}/{traits}/master_table.csv",
            dist_json = "{base}/analyses/{tax_level}/{traits}/tree/rep_dist.json"
    output : pairwise_csv = "{base}/analyses/{tax_level}/{traits}/pairs_table.csv"
    params : script = "workflow/scripts/pairwise.py"
    conda : "../envs/python_std.yaml"
    threads : 1
    shell : """
        python3 {params.script} {input.bmft} {input.dist_json} {output.pairwise_csv} {wildcards.tax_level} {wildcards.traits}
        """
