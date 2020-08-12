library(ggplot2)
library(data.table)
library(hues)
library(vegan)

args = commandArgs(trailingOnly=TRUE)

bmft = fread(args[1])
pairs = fread(args[2])

fig1_file = args[3]
tax_level = args[4]
traits = args[5]

bmft = bmft[mean_genome_size < 15000000]
choice_phylum = names(sort(table(factor(bmft$phylum)), decreasing=TRUE)[1:12])
all_phylum = names(sort(table(factor(bmft$phylum))))
pairs = pairs[sample(nrow(pairs))]
bmft = bmft[sample(nrow(bmft))]


anon_dat = copy(bmft[phylum %in% choice_phylum])
anon_dat[,class := NULL]
anon_dat[,phylum := NULL]

tt = bmft[, .(class)]# %in% choice_phylum, .(class)]
tt[,phylum:=sapply(lapply(strsplit(class,";"), "[", 1:2), paste0, collapse=";")]

class2phylum = tt$phylum
names(class2phylum) = tt$class
class2phylum = class2phylum[!duplicated(names(class2phylum))]
phylum2class = sapply(all_phylum, function(x) names(class2phylum)[class2phylum == x])

base_cols = seq(0,360, 360/(length(all_phylum)+1))[1:length(all_phylum)]
incr = as.numeric(base_cols[2])

base_cols = sample(base_cols)
names(base_cols) =  all_phylum
tt = lapply(all_phylum, function(x){
    n = length(phylum2class[[x]])
    if(n > 1)
    {
      cols = iwanthue(n, hmin = base_cols[x], hmax = base_cols[x]+(incr/2), lmin = 30)
    }
    else
    {
      cols = iwanthue(2, hmin = base_cols[x], hmax = base_cols[x]+(incr/2), lmin = 30)[1]
    }
    names(cols)  = phylum2class[[x]]
    cols
})
class2col = unlist(tt)

setkey(bmft, "clade_name")

annots = fread("analyses/species/ab_initio/annotations/cogcat.csv")

functional_mds = function(blacklist_functions = c('S', '1'), blacklist_clades = NULL, fraction = NA, taxonomy = NA){
  tt = copy(annots[!V1 %in% blacklist_clades])
  if(!is.na(fraction)){
    tt = tt[grepl(fraction,V1)]
  }
  if(!is.na(taxonomy)){
    tt = tt[grepl(taxonomy, bmft[sapply(strsplit(V1,":"), "[" , 1), species])]
  }
  clades = tt$V1

  tt[, V1 := NULL]
  mds = metaMDS(tt)


  data = data.table(mds$points)
  data$clade_name = sapply(strsplit(clades,":"), "[" , 1)
  data = merge(data, bmft, by="clade_name")
  data$facet = data$phylum
  data[!facet %in% choice_phylum, facet := "other"]

  anon_dat = copy(data)
  anon_dat[,class := NULL]
  anon_dat[,facet := NULL]

  ggplot(data[mean_genome_size < 8000000], aes(x=MDS1, y=MDS2, col=class, size=mean_genome_size/1000000))+geom_point(data = anon_dat, col = "grey", fill="grey", alpha=0.3)+geom_point()+scale_color_manual(values = class2col)+guides(col=guide_legend(ncol=1))+theme_minimal()+facet_wrap(~facet)
}









ggplot(bmft[phylum %in% choice_phylum], aes(x=mean_cogs, y=mean_variable, col=order, size = nb_genomes))+
  geom_point()+
#  geom_smooth(data = anon_dat, method = "lm", size = 1, alpha = 0.8, se=F)+
#  facet_wrap(~phylum)+
  guides(col=guide_legend(ncol=2))+
  theme_minimal()+
  scale_y_log10()+scale_x_log10()+
  theme(legend.position = "none")+xlim(500,8000)+ylim(500,10000)

ggplot(bmft[phylum %in% choice_phylum], aes(x=mean_genome_size/1000000, y=mean_GC, col=class, size = sqrt(nb_estgenomes),shape = mOTU_type))+
#  geom_point(data = anon_dat[phylum %in% choice_phylum], col="grey", alpha=0.5)+
  geom_point()+
  facet_wrap(~phylum)+#, scale = "free_x")+
  guides(col=guide_legend(ncol=2))+
  theme_minimal()+#scale_x_log10(lim = c(0.5,8))+theme(legend.position = "none")+
  scale_color_manual(values = class2col)

    ggplot(bmft[phylum %in% choice_phylum], aes(x=mean_GC, y=mean_varfract, col=class, size = sqrt(nb_genomes)))+
      geom_point()+
    #  geom_smooth(data = anon_dat, method = "lm", size = 1, alpha = 0.8, se=F)+
  #    facet_wrap(~phylum, scale = "free")+
      guides(col=guide_legend(ncol=2))+
      theme_minimal()+scale_y_log10()+
      theme(legend.position = "none")+scale_color_manual(values = class2col)


ggplot(bmft[phylum %in% choice_phylum], aes(x=mean_cogs, y=accessory_len/mean_variable/log(nb_estgenomes), col=order, size = nb_genomes))+
  geom_point()+
#  geom_smooth(data = anon_dat, method = "lm", size = 1, alpha = 0.8, se=F)+
  facet_wrap(~phylum)+
  guides(col=guide_legend(ncol=2))+
  theme_minimal()+
  scale_y_log10()+scale_x_log10(lim = c(500,4000))+
  theme(legend.position = "none")

ggplot(bmft[phylum %in% choice_phylum], aes(x=mean_cogs, y=accessory_len/mean_variable/sqrt(nb_genomes), col=class, size = nb_genomes))+
  geom_point()+
#  geom_smooth(data = anon_dat, method = "lm", size = 1, alpha = 0.8, se=F)+
  facet_wrap(~phylum)+
  guides(col=guide_legend(ncol=2))+
  theme_minimal()+scale_x_log10()+
  theme(legend.position = "none")+xlim(500,7000)


ggplot(pairs[phylum %in% choice_phylum & 10^(cogs_ratio) <100], aes(x = phylo_dist , y = 10^(varfract_ratio), col=class))+geom_point(alpha = 0.5)+scale_color_manual(values = class2col)+theme_minimal()+scale_y_continuous(trans="log2")+scale_x_log10()+facet_wrap(~phylum)
