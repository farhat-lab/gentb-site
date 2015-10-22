library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(countrycode)

##----- Phenotype information
# From the Dataverse:
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/GYHIB2&version=DRAFT

load("Strain_Host_Phenotype.RData")
dv_phenotype <- tbl_dt(data6)
rm(data6)

##----- Clean country names

# sort(unique(dv_phenotype$Country))

dv_phenotype$Country <- as.character(dv_phenotype$Country)
dv_phenotype$Country[dv_phenotype$Country == "unknown"] <- NA
dv_phenotype$Country[dv_phenotype$Country == "China -Tibet"] <- "Tibet"
dv_phenotype$Country[dv_phenotype$Country == "montenegro"] <- "Montenegro"
dv_phenotype$Country[dv_phenotype$Country == "Philipines"] <- "Philippines"
dv_phenotype$Country[dv_phenotype$Country == "South Afica"] <- "South Africa"
dv_phenotype$Country[dv_phenotype$Country == "South Korea N"] <- "Korea"
dv_phenotype$Country[dv_phenotype$Country == "Netherlands Antilles"] <- "Netherlands Antilles"
dv_phenotype$Country[dv_phenotype$Country == "Ivory Coast"] <- "Cote d'Ivoire"
dv_phenotype$Country[dv_phenotype$Country == "Guinea Eq."] <- "Equatorial Guinea"
dv_phenotype$Country[dv_phenotype$Country == "Guinea-Conakry"] <- "Guinea"
dv_phenotype$Country[dv_phenotype$Country == "PR"] <- "Puerto Rico"
dv_phenotype$Country[dv_phenotype$Country == "Comoro Islands"] <- "Comoros"
dv_phenotype$Country[dv_phenotype$Country == "Dominican Republic"] <- "Dominican Republic"
dv_phenotype$Country[dv_phenotype$Country == "RD Congo"] <- "Congo, the Democratic Republic of the"
dv_phenotype$Country[dv_phenotype$Country == "DR Congo"] <- "Congo, the Democratic Republic of the"
dv_phenotype$Country[dv_phenotype$Country == "Marocco"] <- "Morocco"

# sort(unique(dv_phenotype$Country))

##----- Clean family names

# sort(unique(dv_phenotype$Spfamily_parentstrain))

dv_phenotype$Spfamily_parentstrain[dv_phenotype$Spfamily_parentstrain == "beijing"] <- "Beijing"
dv_phenotype$Spfamily_parentstrain[dv_phenotype$Spfamily_parentstrain == "BEIJING"] <- "Beijing"
dv_phenotype$Spfamily_parentstrain[dv_phenotype$Spfamily_parentstrain == "Beijing-1-atypical"] <- "Beijing"
dv_phenotype$Spfamily_parentstrain[dv_phenotype$Spfamily_parentstrain == "beijing-6-typical"] <- "Beijing"
# dv_phenotype$Spfamily_parentstrain[dv_phenotype$Spfamily_parentstrain == "CDC 1551"] <- "CDC"
# dv_phenotype$Spfamily_parentstrain[dv_phenotype$Spfamily_parentstrain == "CDC-609"] <- "CDC"
# dv_phenotype$Spfamily_parentstrain[dv_phenotype$Spfamily_parentstrain == "CDC-610"] <- "CDC"

# sort(unique(dv_phenotype$Spfamily_parentstrain))

##----- Create country codes

dv_phenotype <- dv_phenotype[!is.na(Country)]
dv_phenotype[, Country_code := countrycode(Country, "country.name", "iso3c")]

##----- Country centroids with country code information

country_geoinformation <- read.csv("Country_List_ISO_3166_Codes_Latitude_Longitude.csv")
country_geoinformation <- country_geoinformation[c("Alpha.3.code", "Latitude..average.", "Longitude..average.")]
names(country_geoinformation) <- c("Country_code", "lat", "lon")

##----- Genotype information

# From the Dataverse:
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/GYHIB2&version=DRAFT
load("matrix.RData")

dv_genotype <- tbl_dt(matrix)
rm(matrix)

dv_genotype_vstack_sparse <- dv_genotype %>%
  gather(strain, value)
setnames(dv_genotype_vstack_sparse, "strain.1", "snpname")

dv_genotype_vstack <- subset(dv_genotype_vstack_sparse, value == 1)
dv_genotype_vstack[, value := NULL]
dv_genotype_vstack[, snpname := as.character(snpname)]

d_split_name <- tbl_df(as.data.frame(str_split_fixed(dv_genotype_vstack$snpname, "_", 6)))
colnames(d_split_name) <- c("sub.ins.del",
                            "coding.silent",
                            "genomic.coordinate",
                            "genetic.change",
                            "amino.acid.change",
                            "region")

dt_split_name <- data.table(d_split_name)
dt_split_name[amino.acid.change == "promoter", region := paste0("promoter_", region)]
dt_split_name[amino.acid.change == "promoter", amino.acid.change := "NA"]
dt_split_name[amino.acid.change == "inter", region := paste0("inter_", region)]
dt_split_name[amino.acid.change == "inter", amino.acid.change := "NA"]

dv_genotype_vstack <- cbind(dv_genotype_vstack, dt_split_name)

##----- Merge datasets

## All

dv_phenotype_with_country <- merge(dv_phenotype, country_geoinformation, by = "Country_code")

setkey(dv_phenotype_with_country, Strain)
setkey(dv_genotype_vstack, strain)

dv_genotype_phenotype_with_country <- dv_genotype_vstack[dv_phenotype_with_country]
write.csv(dv_genotype_phenotype_with_country, "dv_genotype_phenotype_with_country.csv")

## By strain

dv_genotype_phenotype_with_country_by_strain <- dplyr::select(tbl_df(unique(dv_genotype_phenotype_with_country, by = "strain")),
                                                              strain, ResistanceType, ends_with("DST", ignore.case = FALSE),
                                                              Country, Country_code, lat, lon,
                                                              Spfamily_parentstrain, RFLPfamily)
write.csv(dv_genotype_phenotype_with_country_by_strain, "dv_genotype_phenotype_with_country_by_strain.csv")

#----- Unique by patient

# m <- data.table(m)
# setkey(m, ID)
# setkeyv(m, c("code", "ID")) # sort by country and then by patient
# mu <- unique(m)
# 
# write.csv(mu[, list(ID, drtype,
#                     rinh, reth, rrif, rrfb, remb,
#                     rpza, rstr, ramk, rcap,
#                     rkan, roflx, rcip, rlevo,
#                     rpas, rcys, rtha, rpro,
#                     rclof, rmoxi, rclar, rgati,              
#                     ramoxclav, rlin,
#                     spfamily_parentstrain,
#                     code, country, lat, lon)], "merged_unique_patient.csv")

#----- TO REMOVE?

# library(dplyr)
# pu <- tbl_dv_phenotype(read.csv("merged_unique_patient.csv"))
# ppu <- group_by(pu, code)
# 
# pppu <- group_by(ppu, drtype)
# summarize(pppu, cnt = n())
# 
# ppur <- select(ppu, starts_with("r"))
# summarize(ppur, cnt = n())
# 
# pu %>%
#   group_by(code) %>%
#   summarize(
#     country = country,
#     strain = n(),
#     lat = lat,
#     lon = lon
#     ) %>%
#   distinct()
# 
# ppur <- summarize(ppu,
#           country = country,
#           strain = n(),
#           lat = lat,
#           lon = lon
#           )

#----- Create Shiny dataset

# d.strain <- cbind(table(mu$code))
# colnames(d.strain) <- "strain"
# 
# d.drtype <- table(mu$code, mu$drtype)
# colnames(d.drtype) <- paste0("drtype.", colnames(d.drtype))
# 
# antibio <- c("rinh", "reth", "rrif", "rrfb",                 
#              "remb", "rpza", "rstr", "ramk",
#              "rcap", "rkan", "roflx", "rcip", 
#              "rlevo", "rpas", "rcys", "rtha",
#              "rpro", "rclof", "rmoxi", "rclar",
#              "rgati", "ramoxclav", "rlin")
# 
# d.antibio <- NULL
# 
# for (a in antibio) {
#   d.temp <- table(mu$code, mu[[a]])
#   colnames(d.temp) <- paste0(a, ".", colnames(d.temp))
#   d.antibio <- cbind(d.antibio, d.temp)
# }
# 
# d.spfamily_parentstrain <- table(mu$code, mu$spfamily_parentstrain)
# colnames(d.spfamily_parentstrain) <- paste0("family.", colnames(d.spfamily_parentstrain))
# 
# m.shiny <- cbind(d.strain,
#                  d.drtype,
#                  d.antibio,
#                  d.spfamily_parentstrain)
# 
# mm.shiny <- merge(m.shiny, country_geoinformation)
# 
# write.csv(m.shiny, "m_shiny.csv")
# 
# 

