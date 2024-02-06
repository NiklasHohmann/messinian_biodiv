messinian_db = read.csv(file = "data/messinianDB.csv")

time_bins = unique(messinian_db$Age)[c(3,1,2)]
regions = unique(messinian_db$region.new)
group_names = unique(messinian_db$group.name)
spec = c("occ_all", "occ_eMed", "occ_wMed", "occ_PoA", "species", "genera", "families", "loc_all", "loc_eMed", "loc_wMed", "loc_PoA")

table = matrix(data = NA,
               nrow = length(group_names) + 1,
               ncol = length(spec),
               dimnames = list("group" = c(group_names, "total"),
                               "spec" = spec))
table["benthic_foraminifera","loc_PoA"]
