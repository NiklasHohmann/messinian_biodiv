cat("Running analysis. This might take a few minutes, please wait...\n")

#### load data ####
cat("Loading data\n")
messinian_db <- read.csv(file = "data/messinianDB.csv")

#### Clean data ####
remove_occ = grepl("Considered reworked", messinian_db$Notes, useBytes = TRUE) | grepl("Collected from deposits known as \"Livelli ad Aturia", messinian_db$Notes, useBytes = TRUE)
messinian_db = messinian_db[!remove_occ,]

#### clean coral type data ####
corals_data = read.csv(file = "data/coral genera FB.csv")
corals_data$type = rep(NA, length(corals_data$Genus))
corals_data$type = replace(corals_data$type, corals_data$Zooxanthellate == "x", "z")
corals_data$type = replace(corals_data$type, corals_data$Azooxanthellate == "x", "a")
corals_data$type = replace(corals_data$type, corals_data$notes == "in reef setting", "z")
corals_data$type = replace(corals_data$type, corals_data$notes == "not in reef setting", "a")

# introduce new tax group classification where a and z coral are analyzed separately
messinian_db$group.name_sc = messinian_db$group.name  # sc = split corals
for (i in seq_along(messinian_db$ID)){
  if (messinian_db$group.name[i] == "corals" & messinian_db$Genus.name[i] %in% corals_data$Genus){
    type = corals_data$type[corals_data$Genus == messinian_db$Genus.name[i]]
    messinian_db$group.name_sc[i] = paste(messinian_db$group.name[i], type, sep = "_")
  }
}

#### Set seed ####
set.seed(1)

#### Replace invalid tax names with NA ####
messinian_db$Species.name = replace(messinian_db$Species.name, messinian_db$Species.name == "sp.", NA)
messinian_db$Genus.name = replace(messinian_db$Genus.name, messinian_db$Genus.name == "indet.", NA)
messinian_db = messinian_db[ !messinian_db$Family == "indet.",  ]



#### Constants ####
timebins <- unique(messinian_db$Age)[c(3, 1, 2)] # sorted from old to young
regions <- unique(messinian_db$region.new)
group.names <- unique(messinian_db$group.name)
group.names.ext = c(group.names, "all groups")
group.names_sc  = unique(messinian_db$group.name_sc[messinian_db$group.name_sc != "corals"]) # with corals split
group.names.ext_sc = c(group.names_sc, "all groups") 
regions.ext = c(regions, "whole basin")
eco_index_names = c("soerensen", "simpson", "nestedness")
timebin_comp = c("T vs. M", "M vs. Z", "T vs. Z")
noOfRep = 10000 # number of repetitions for subsampling
taxonomic_levels = c("species", "genus")

#### Load helper functions ####
cat("Loading helper functions\n")
source("code/helper_functions.R")

#### Taxonomic richness through time (whole basin, all groups) ####
cat("Determining taxonomic richness\n")
for (tax_level in taxonomic_levels){
  dir.create(paste0("figs/sr_through_time/",tax_level), showWarnings = FALSE)
  for (group in group.names.ext){
    # extract species names
    Tor = get_from_db(group = group, basin = "whole basin", timeslice = "Tortonian" , taxLevel = tax_level)
    Mes = get_from_db(group = group, basin = "whole basin", timeslice = "pre-evaporitic Messinian" , taxLevel = tax_level)
    Zan = get_from_db(group = group, basin = "whole basin", timeslice = "Zanclean", taxLevel = tax_level)
    # define subsampling size (80 % of smallest sample)
    subsampleTo = ceiling(0.8 * (min(c(length(Tor), length(Mes), length(Zan)))))
    # subsample noOfRep times
    Tor_sr = rarefyTaxRichness(mySample = Tor, subsampleTo = subsampleTo, noOfRep = noOfRep)
    Mes_sr = rarefyTaxRichness(mySample = Mes, subsampleTo = subsampleTo, noOfRep = noOfRep)
    Zan_sr = rarefyTaxRichness(mySample = Zan, subsampleTo = subsampleTo, noOfRep = noOfRep)
    # make figure
    file_name = paste0("figs/sr_through_time/", tax_level, "/sr_through_time_whole_basin_",group, "_", tax_level, "level.pdf")
    main = paste0("Species Richness ", group, " whole basin, ", tax_level, " level")
    ylim = c(0, max(c(Tor_sr, Mes_sr, Zan_sr)))
    ylab = paste0("Species richness \n subsampled to ",  subsampleTo, " Occurrences")
    pdf(file = file_name)
    boxplot(list( "Tortonian" = Tor_sr,
                  "Messinian" = Mes_sr,
                  "Zanclean" =  Zan_sr),
            ylim = ylim,
            ylab = ylab,
            main = main)
    dev.off()
  }
}

for (tax_level in c("genus")){
  for (group in c("corals_a", "corals_z")){
    # extract species names
    Tor = get_from_db(group = group, basin = "whole basin", timeslice = "Tortonian" , tax_groups = "split corals", taxLevel = tax_level)
    Mes = get_from_db(group = group, basin = "whole basin", timeslice = "pre-evaporitic Messinian" , tax_groups = "split corals", taxLevel =  tax_level)
    Zan = get_from_db(group = group, basin = "whole basin", timeslice = "Zanclean" , tax_groups = "split corals", taxLevel = tax_level)
    # define subsampling size (80 % of smallest sample)
    subsampleTo = ceiling(0.8 * (min(c(length(Tor), length(Mes), length(Zan)))))
    # subsample noOfRep times
    Tor_sr = rarefyTaxRichness(mySample = Tor, subsampleTo = subsampleTo, noOfRep = noOfRep)
    Mes_sr = rarefyTaxRichness(mySample = Mes, subsampleTo = subsampleTo, noOfRep = noOfRep)
    Zan_sr = rarefyTaxRichness(mySample = Zan, subsampleTo = subsampleTo, noOfRep = noOfRep)
    # make figure
    file_name = paste0("figs/sr_through_time/", tax_level, "/sr_through_time_whole_basin_",group, "_", tax_level, "_level.pdf")
    main = paste0("Species Richness ", group, " whole basin, ",  tax_level, " level")
    ylim = c(0, max(c(Tor_sr, Mes_sr, Zan_sr)))
    ylab = paste0("taxonomic richness \n subsampled to ",  subsampleTo, " Occurrences ")
    pdf(file = file_name)
    boxplot(list( "Tortonian" = Tor_sr,
                  "Messinian" = Mes_sr,
                  "Zanclean" =  Zan_sr),
            ylim = ylim,
            ylab = ylab,
            main = main)
    dev.off()
  }
}


#### Ecological indices through time (whole basin, all groups) ####
cat("Determining ecological indices\n")
for (tax_level in taxonomic_levels){
  for (group in group.names.ext){
    dir.create(paste0("figs/eco_timeslice_comp/", tax_level), showWarnings = FALSE)
    # extract species names
    Tor = get_from_db(group = group, basin = "whole basin", timeslice = "Tortonian", taxLevel = tax_level )
    Mes = get_from_db(group = group, basin = "whole basin", timeslice = "pre-evaporitic Messinian" , taxLevel = tax_level)
    Zan = get_from_db(group = group, basin = "whole basin", timeslice = "Zanclean" , taxLevel = tax_level)
    # define subsampling size (80 % of smallest sample)
    subsampleTo = ceiling(0.8 * (min(c(length(Tor), length(Mes), length(Zan)))))
    # subsample noOfRep times
    TM = rarefyEcoIndexes(Tor, Mes, subsampleTo, noOfRep)
    MZ = rarefyEcoIndexes(Mes, Zan, subsampleTo, noOfRep)
    TZ = rarefyEcoIndexes(Tor, Zan, subsampleTo, noOfRep)
    # make plots
    for (ind in eco_index_names){
      file_name = paste0("figs/eco_timeslice_comp/",tax_level, "/", group, "_", ind, "_", tax_level, "level_whole_basin.pdf")
      pdf(file = file_name)
      main = paste0(ind, " for " , group, " whole basin ", tax_level, "level")
      ylab = paste0(ind, "\n subsampled to ",  subsampleTo, " Occurrences")
      boxplot(list("T vs. M" = TM[[ind]],
                   "M vs. Z" = MZ[[ind]],
                   "T vs. Z" = TZ[[ind]]),
              ylim = c(0,1),
              main = main,
              ylab = ylab,
              mar = c(5,5,1,1))
      dev.off()
    }
  }
}

for (tax_level in c("genus")){
  for (group in c("corals_a", "corals_z")){
    # extract species names
    Tor = get_from_db(group = group, basin = "whole basin", timeslice = "Tortonian", tax_groups = "split corals" , taxLevel = tax_level)
    Mes = get_from_db(group = group, basin = "whole basin", timeslice = "pre-evaporitic Messinian" , tax_groups = "split corals", taxLevel = tax_level)
    Zan = get_from_db(group = group, basin = "whole basin", timeslice = "Zanclean" , tax_groups = "split corals", taxLevel = tax_level)
    # define subsampling size (80 % of smallest sample)
    subsampleTo = ceiling(0.8 * (min(c(length(Tor), length(Mes), length(Zan)))))
    # subsample noOfRep times
    TM = rarefyEcoIndexes(Tor, Mes, subsampleTo, noOfRep)
    MZ = rarefyEcoIndexes(Mes, Zan, subsampleTo, noOfRep)
    TZ = rarefyEcoIndexes(Tor, Zan, subsampleTo, noOfRep)
    # make plots
    for (ind in eco_index_names){
      file_name = paste0("figs/eco_timeslice_comp/",tax_level, "/", group, "_", ind,"_", tax_level, "_level_whole_basin.pdf")
      pdf(file = file_name)
      main = paste0(ind, " for " , group, " whole basin")
      ylab = paste0(ind, "\n subsampled to ",  subsampleTo, " Occurrences", tax_level, "level")
      boxplot(list("T vs. M" = TM[[ind]],
                   "M vs. Z" = MZ[[ind]],
                   "T vs. Z" = TZ[[ind]]),
              ylim = c(0,1),
              main = main,
              ylab = ylab,
              mar = c(5,5,1,1))
      dev.off()
    }
  }
}

#### Taxonomic richness through time for all groups in all regions ####
cat("Determining regional species richness\n")
sr_median = array(data = NA, # median species richness
                       dim = c(length(taxonomic_levels), length(group.names.ext), length(timebins), length(regions.ext)),
                       dimnames = list("tax_level" = taxonomic_levels,
                                        "group" = group.names.ext,
                                       "timebin" = timebins,
                                       "region" = regions.ext))
for(tax_level in taxonomic_levels){
  dir.create(paste0("figs/sr_through_time_regional/", tax_level), showWarnings = FALSE)
  for (group in group.names.ext){
    for (reg in regions.ext){
      Tor = get_from_db(group = group, basin = reg, timeslice = "Tortonian", taxLevel = tax_level )
      Mes = get_from_db(group = group, basin = reg, timeslice = "pre-evaporitic Messinian" , taxLevel = tax_level)
      Zan = get_from_db(group = group, basin = reg, timeslice = "Zanclean" , taxLevel = tax_level)
      subsampleTo = ceiling(0.8 * (min(c(length(Tor), length(Mes), length(Zan)))))
      Tor_sr = rarefyTaxRichness(mySample = Tor, subsampleTo = subsampleTo, noOfRep = noOfRep)
      Mes_sr = rarefyTaxRichness(mySample = Mes, subsampleTo = subsampleTo, noOfRep = noOfRep)
      Zan_sr = rarefyTaxRichness(mySample = Zan, subsampleTo = subsampleTo, noOfRep = noOfRep)
      sr_median[tax_level, group, "Tortonian", reg] = median(Tor_sr)
      sr_median[tax_level, group, "pre-evaporitic Messinian", reg] = median(Mes_sr)
      sr_median[tax_level, group, "Zanclean", reg] = median(Zan_sr)
      
      file_name = paste0("figs/sr_through_time_regional/", tax_level, "/sr_through_time_", group, "_",  reg , "_", tax_level, "_level.pdf")
      ylab = paste0("Species richness \n subsampled to ",  subsampleTo, " Occurrences")
      main = paste0("Species Richness ", group, " ", reg, " ", tax_level, " level")
      ylim = c(0, max(c(Tor_sr, Mes_sr, Zan_sr)))
      pdf(file = file_name)
      boxplot(list( "Tortonian" = Tor_sr,
                    "Messinian" = Mes_sr,
                    "Zanclean" =  Zan_sr),
              ylim = ylim,
              ylab = ylab,
              main = main)
      dev.off()
    }
  }
}

#### Ecological indices through time per region ( all groups) ####

eco_ind_median = array(data = NA,
                         dim = c(length(taxonomic_levels), length(group.names.ext), length(timebin_comp), length(regions.ext), length(eco_index_names)),
                         dimnames = list("tax_level" = taxonomic_levels,
                                         "group" = group.names.ext,
                                         "timebins" = timebin_comp,
                                         "regions" = regions.ext,
                                         "index" = eco_index_names))

cat("Determining regional ecological indices\n")
for (tax_level in taxonomic_levels){
  dir.create(paste0("figs/eco_timeslice_comp_regional/", tax_level), showWarnings = FALSE)
  for (group in c(group.names.ext)){
    for (reg in regions.ext){
      Tor = get_from_db(group = group, basin = reg, timeslice = "Tortonian", taxLevel = tax_level )
      Mes = get_from_db(group = group, basin = reg, timeslice = "pre-evaporitic Messinian" , taxLevel = tax_level)
      Zan = get_from_db(group = group, basin = reg, timeslice = "Zanclean", taxLevel = tax_level )
      subsampleTo = ceiling(0.8 * (min(c(length(Tor), length(Mes), length(Zan)))))
      TM = rarefyEcoIndexes(Tor, Mes, subsampleTo, noOfRep)
      MZ = rarefyEcoIndexes(Mes, Zan, subsampleTo, noOfRep)
      TZ = rarefyEcoIndexes(Tor, Zan, subsampleTo, noOfRep)
      for (ind in eco_index_names){
        eco_ind_median[ tax_level, group, "T vs. M", reg, ind] = 100 * median(TM[[ind]])
        eco_ind_median[tax_level, group, "M vs. Z", reg, ind] = 100 * median(MZ[[ind]])
        eco_ind_median[tax_level, group, "T vs. Z", reg, ind] = 100 * median(TZ[[ind]])
        file_name = paste0("figs/eco_timeslice_comp_regional/",tax_level, "/", group, "_", reg ,"_", ind, "_", tax_level, ".pdf")
        pdf(file = file_name)
        main = paste0(ind, " for " , group, " ", reg, " ", tax_level)
        ylab = paste0(ind, "\n subsampled to ",  subsampleTo, " Occurrences")
        boxplot(list("T vs. M" = TM[[ind]],
                     "M vs. Z" = MZ[[ind]],
                     "T vs. Z" = TZ[[ind]]),
                ylim = c(0,1),
                main = main,
                ylab = ylab,
                mar = c(5,5,1,1))
        dev.off()
      }
    }
  }
}

#### Taxonomic richness per group and region, rarefied to the same sample size ####
cat("Determining regional species richness, rarefied to same sample size\n")
sr_median_comparable = array(data = NA, # median species richness
                            dim = c(length(taxonomic_levels), length(group.names.ext), length(timebins), length(regions.ext)),
                            dimnames = list("tax_level" = taxonomic_levels,
                                            "group" = group.names.ext,
                                            "timebin" = timebins,
                                            "region" = regions.ext))
for (tax_level in taxonomic_levels){
  dir.create(paste0("figs/sr_through_time_regional_comparable/", tax_level), showWarnings = FALSE)
  for (group in group.names.ext){
    # determine subsampling size
    subsampleTo = Inf
    for (ti in timebins){
      wMed = get_from_db(group, "Western Mediterranean", ti, taxLevel = tax_level)
      eMed = get_from_db(group, "Eastern Mediterranean", ti, taxLevel = tax_level)
      PPNA = get_from_db(group, "Po Plain-Northern Adriatic", ti, taxLevel = tax_level)
      subsampleTo = min(subsampleTo, ceiling(0.8 * (min(c(length(wMed), length(eMed), length(PPNA))))))
    }
    for (reg in regions.ext){
      Tor = get_from_db(group = group, basin = reg, timeslice = "Tortonian" , taxLevel = tax_level)
      Mes = get_from_db(group = group, basin = reg, timeslice = "pre-evaporitic Messinian" , taxLevel = tax_level)
      Zan = get_from_db(group = group, basin = reg, timeslice = "Zanclean" , taxLevel = tax_level)
      Tor_sr = rarefyTaxRichness(mySample = Tor, subsampleTo = subsampleTo, noOfRep = noOfRep)
      Mes_sr = rarefyTaxRichness(mySample = Mes, subsampleTo = subsampleTo, noOfRep = noOfRep)
      Zan_sr = rarefyTaxRichness(mySample = Zan, subsampleTo = subsampleTo, noOfRep = noOfRep)
      sr_median_comparable[tax_level, group, "Tortonian", reg] = median(Tor_sr)
      sr_median_comparable[tax_level, group, "pre-evaporitic Messinian", reg] = median(Mes_sr)
      sr_median_comparable[tax_level, group, "Zanclean", reg] = median(Zan_sr)
      
      file_name = paste0("figs/sr_through_time_regional_comparable/", tax_level, "/sr_through_time_comp", group, "_",  reg , "_", tax_level, ".pdf")
      ylab = paste0("Species richness \n subsampled to ",  subsampleTo, " Occurrences")
      main = paste0("Species Richness ", group, " ", reg, " ", tax_level)
      ylim = c(0, max(c(Tor_sr, Mes_sr, Zan_sr)))
      pdf(file = file_name)
      boxplot(list( "Tortonian" = Tor_sr,
                    "Messinian" = Mes_sr,
                    "Zanclean" =  Zan_sr),
              ylim = ylim,
              ylab = ylab,
              main = main)
      dev.off()
    }
  }
}


#### Changes in sr & eco indexes ####
sr_change = array(data = NA,
                  dim = c(length(taxonomic_levels), length(group.names.ext), 3, length(regions.ext)),
                  dimnames = list("tax_level" = taxonomic_levels, 
                                  "group" = group.names.ext,
                                   "time_slices" = c("TvsM","MvsZ","TvsZ"),
                                  "region" = regions.ext))
for (tax_level in taxonomic_levels){
  for (group in group.names.ext){
    for (reg in regions.ext){
      sr_change[tax_level, group, "TvsM", reg] = 100 * (1 - sr_median[tax_level, group, "Tortonian",reg] / sr_median[tax_level, group,"pre-evaporitic Messinian",reg])
      sr_change[tax_level, group, "MvsZ", reg] = 100 * (1 - sr_median[tax_level, group, "pre-evaporitic Messinian",reg] / sr_median[tax_level, group,"Zanclean",reg])
      sr_change[tax_level, group, "TvsZ", reg] = 100 * (1 - sr_median[tax_level, group, "Tortonian",reg] / sr_median[tax_level, group,"Zanclean",reg])
    }
  }
}

cat("Done. Outputs are in the folder \"figs\" and the variables \"sr_change\" and \"eco_ind_median\".\n")

## assessing sr changes in PERCENT
# first entry: taxonomic level out of taxonomic_levels
# second entry: name out of group.names.ext.
# third entry: "TvsM","MvsZ", or "TvsZ": time slice to compare
# fourth entry: name out of regions.ext: region to examine
# returns change in median species richness
#sr_change["species", "all groups", "TvsM", "whole basin"] 

## extracting median ecological indices in PERCENT
# first entry: taxonomic level out of taxonomic_levels
# second entry: group name, element of group.names.ext
# third entry: time bins to compare, "T vs. M", "M vs. Z" or "T vs. Z"
# fourth entry: region of interest, element of regions.ext
# fifth entry: ecological index, element of eco_index_names
# returns median ecological index for the specified case.
#eco_ind_median["species",  "fish", "T vs. M", "whole basin", "simpson"]


# sr_change["species",,"TvsM", "whole basin"]
# sr_change["species",,"MvsZ", "whole basin"]
# eco_ind_median["species",, "T vs. M", "whole basin", "soerensen"]
# eco_ind_median["species",, "M vs. Z", "whole basin", "soerensen"]
# eco_ind_median["species",, "T vs. Z", "whole basin", "soerensen"]
# eco_ind_median["species",, "T vs. Z", "whole basin", "simpson"]
# eco_ind_median["species",, "M vs. Z", "whole basin", "simpson"]

#### Test for statements in ms ####

# genus richness of bivalves shows a smaller decrease from the Tortonian to the pre-evaporitic Messinian
group.name = "bivalves"
t1 = timebins[1]
t2 = timebins[2]

biv_tor_gen = get_from_db(group.name, "whole basin", t1, "genus")
biv_mes_gen = get_from_db(group.name, "whole basin", t2, "genus")
biv_tor_spe = get_from_db(group.name, "whole basin", t1)
biv_mes_spe = get_from_db(group.name, "whole basin", t2)

subsampleTo = min(ceiling(0.8 * c(length(biv_tor_gen), length(biv_mes_gen), length(biv_tor_spe), length(biv_mes_spe))))
biv_gen_gradient = rarefyTaxGradient(biv_tor_gen, biv_mes_gen, subsampleTo, noOfRep)
biv_spe_gradient = rarefyTaxGradient(biv_tor_spe, biv_mes_spe, subsampleTo, noOfRep)

biv_diff_tax_grad_test = wilcox.test(biv_gen_gradient, biv_spe_gradient, alternative = "less")

# Gastropod genus richness shows marginal increase in the Zanclean (Fig. S1H), while species richness remains the same 
group.name = "gastropods"
t1 = timebins[2]
t2 = timebins[3]

gas_mes_gen = get_from_db(group.name, "whole basin", t1, "genus")
gas_zan_gen = get_from_db(group.name, "whole basin", t2, "genus")
gas_mes_spe = get_from_db(group.name, "whole basin", t1)
gas_zan_spe = get_from_db(group.name, "whole basin", t2)

subsampleTo = min(ceiling(0.8 * c(length(gas_mes_gen), length(gas_mes_gen), length(gas_zan_spe), length(gas_zan_spe))))

gen_grad = rarefyTaxGradient(gas_mes_gen, gas_zan_gen, subsampleTo, noOfRep)
spe_grad = rarefyTaxGradient(gas_mes_spe, gas_zan_spe, subsampleTo, noOfRep)

gas_gen_test_dec = wilcox.test(gen_grad, alternative = "less")
gas_spe_test_unequal = wilcox.test(spe_grad)
quantile(spe_grad)


# For echinoids, genus richness increases (Fig. S1J), but species richness decreases (Fig. 3J) from the Tortonian to the pre-evaporitic Messinian

group.name = "echinoids"
t1 = timebins[1]
t2 = timebins[2]

ech_tor_gen = get_from_db(group.name, "whole basin", t1, "genus")
ech_mes_gen = get_from_db(group.name, "whole basin", t2, "genus")
ech_tor_spe = get_from_db(group.name, "whole basin", t1)
ech_mes_spe = get_from_db(group.name, "whole basin", t2)

subsampleTo = min(ceiling(0.8 * c(length(ech_mes_gen), length(ech_mes_gen), length(ech_tor_spe), length(ech_mes_spe))))

gen_grad = rarefyTaxGradient(ech_tor_gen, ech_mes_gen, subsampleTo, noOfRep)
spe_grad = rarefyTaxGradient(ech_tor_spe, ech_mes_spe, subsampleTo, noOfRep)

ech_gen_inc_test = wilcox.test(gen_grad , alternative = "less")
ech_spe_dec_test = wilcox.test(spe_grad, alternative = "greater")

# Species richness decreased from the Tortonian to the pre-evaporitic Messinian for calcareous nannoplankton, dinocysts, , planktic foraminifera , corals , bivalves, echinoids, and bony fishes 
sel_groups = c("nanoplankton", "dinocysts", "planktic_foraminifera", "corals", "bivalves", "echinoids", "fish")
t1 = timebins[1]
t2 = timebins[2]
testreslist1  = list()
for (group.name in sel_groups){
  sp_r_tor = get_from_db(group = group.name, basin = "whole basin", timeslice = t1)
  sp_r_mes = get_from_db(group = group.name, basin = "whole basin", timeslice = t2)
  sp_r_zan = get_from_db(group = group.name, basin = "whole basin", timeslice = timebins[3])
  subsampleTo = ceiling(0.8 * (min(c(length(sp_r_tor), length(sp_r_mes), length(sp_r_zan)))))
  sp_grad = rarefyTaxGradient(sp_r_tor, sp_r_mes, subsampleTo, noOfRep)
  testreslist1[[group.name]] = wilcox.test(sp_grad, alternative = "greater")
}
p_vals = sapply(testreslist1, function(x) x$p.value)

# ostracod species richness remains at the same level

group.name = "ostracods"

t1 = timebins[1]
t2 = timebins[2]
  sp_r_tor = get_from_db(group = group.name, basin = "whole basin", timeslice = t1)
  sp_r_mes = get_from_db(group = group.name, basin = "whole basin", timeslice = t2)
  sp_r_zan = get_from_db(group = group.name, basin = "whole basin", timeslice = timebins[3])
  subsampleTo = ceiling(0.8 * (min(c(length(sp_r_tor), length(sp_r_mes), length(sp_r_zan)))))
  sp_grad = rarefyTaxGradient(sp_r_tor, sp_r_mes, subsampleTo, noOfRep)
wilcox.test(sp_grad)

# species richness increased from the Tortonian to the Messinian for benthic foraminifera, gastropods, and bryozoans 

sel_groups = c("benthic_foraminifera", "gastropods", "bryozoans")
t1 = timebins[1]
t2 = timebins[2]
testreslist2  = list()
for (group.name in sel_groups){
  sp_r_tor = get_from_db(group = group.name, basin = "whole basin", timeslice = t1)
  sp_r_mes = get_from_db(group = group.name, basin = "whole basin", timeslice = t2)
  sp_r_zan = get_from_db(group = group.name, basin = "whole basin", timeslice = timebins[3])
  subsampleTo = ceiling(0.8 * (min(c(length(sp_r_tor), length(sp_r_mes), length(sp_r_zan)))))
  sp_grad = rarefyTaxGradient(sp_r_tor, sp_r_mes, subsampleTo, noOfRep)
  testreslist2[[group.name]] = wilcox.test(sp_grad, alternative = "less")
}
p_vals = sapply(testreslist2, function(x) x$p.value)

# From the Messinian to the Zanclean, species richness remained approximately the same for calcareous nannoplankton planktic and benthic foraminifera, gastropods, and bryozoans

sel_groups = c("nanoplankton", "benthic_foraminifera", "gastropods", "bryozoans"  )
t1 = timebins[2]
t2 = timebins[3]
testreslist3 = list()

for (group.name in sel_groups){
  sp_r_tor = get_from_db(group = group.name, basin = "whole basin", timeslice = timebins[1])
  sp_r_mes = get_from_db(group = group.name, basin = "whole basin", timeslice = t1)
  sp_r_zan = get_from_db(group = group.name, basin = "whole basin", timeslice = t2)
  subsampleTo = ceiling(0.8 * (min(c(length(sp_r_tor), length(sp_r_mes), length(sp_r_zan)))))
  sp_grad = rarefyTaxGradient(sp_r_mes, sp_r_zan, subsampleTo, noOfRep)
  testreslist3[[group.name]] = wilcox.test(sp_grad)
  hist(sp_grad, main = paste( group.name, t1, t2))
}
p_vals = sapply(testreslist3, function(x) x$p.value)

# From the Messinian to the Zanclean, there is a decrease in richness for dinocysts, ostracods, and echinoids. 

sel_groups = c("dinocysts", "ostracods"  ,  "echinoids")
t1 = timebins[2]
t2 = timebins[3]
testreslist4 = list()

for (group.name in sel_groups){
  sp_r_tor = get_from_db(group = group.name, basin = "whole basin", timeslice = timebins[1])
  sp_r_mes = get_from_db(group = group.name, basin = "whole basin", timeslice = t1)
  sp_r_zan = get_from_db(group = group.name, basin = "whole basin", timeslice = t2)
  subsampleTo = ceiling(0.8 * (min(c(length(sp_r_tor), length(sp_r_mes), length(sp_r_zan)))))
  sp_grad = rarefyTaxGradient(sp_r_mes, sp_r_zan, subsampleTo, noOfRep)
  testreslist4[[group.name]] = wilcox.test(sp_grad, alternative = "greater")
  hist(sp_grad, main = paste( group.name, t1, t2))
}
p_vals = sapply(testreslist4, function(x) x$p.value)

# Species richness increased from the Messinian to the Zanclean for corals, bivalves, and bony fishes . 

sel_groups = c("corals", "bivalves"  ,  "fish")
t1 = timebins[2]
t2 = timebins[3]
testreslist5 = list()

for (group.name in sel_groups){
  sp_r_tor = get_from_db(group = group.name, basin = "whole basin", timeslice = timebins[1])
  sp_r_mes = get_from_db(group = group.name, basin = "whole basin", timeslice = t1)
  sp_r_zan = get_from_db(group = group.name, basin = "whole basin", timeslice = t2)
  subsampleTo = ceiling(0.8 * (min(c(length(sp_r_tor), length(sp_r_mes), length(sp_r_zan)))))
  sp_grad = rarefyTaxGradient(sp_r_mes, sp_r_zan, subsampleTo, noOfRep)
  testreslist5[[group.name]] = wilcox.test(sp_grad, alternative = "less")
  hist(sp_grad, main = paste( group.name, t1, t2))
}
p_vals = sapply(testreslist5, function(x) x$p.value)

