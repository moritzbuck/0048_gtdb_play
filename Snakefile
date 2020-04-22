rule proteom_preclustering:
    input : "{taxononomy}/{taxon}.gids":
    output : representatives = "{taxononomy}/{taxon}.preclustering.faa",
             preclusters = "{taxononomy}/{taxon}.preclustering.faa",
    params : script = "scripts/proteom_clustering.py", 
