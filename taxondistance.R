library(chronosphere)
library(divDyn)
library(geosphere)


# Get the data from pbdb
# updated may 17, 2026
dat <- chronosphere::fetch("pbdb", ser="occs4", ver="20260517")

# Filter records not identified at least to genus
dat <- dat[dat$accepted_rank %in% c("genus", "subgenus", "species", "subspecies"),]

# Omit non informative genus entries
dat <- dat[dat$genus != "",]

# Filter phyla
marineNoPlant <- c("",
                   "Agmata",
                   "Annelida",
                   "Bilateralomorpha",
                   "Brachiopoda",
                   "Bryozoa",
                   "Calcispongea",
                   "Chaetognatha",
                   "Cnidaria",
                   "Ctenophora",
                   "Echinodermata",
                   "Entoprocta",
                   "Foraminifera",
                   "Hemichordata",
                   "Hyolitha",
                   "Mollusca",
                   "Nematoda",
                   "Nematomorpha",
                   "Nemertina",
                   "Onychophora",
                   "Petalonamae",
                   "Phoronida",
                   "Platyhelminthes",
                   "Porifera",
                   "Problematica",
                   "Rhizopodea",
                   "Rotifera",
                   "Sarcomastigophora",
                   "Sipuncula",
                   "Uncertain",
                   "Vetulicolia",
                   ""
)

# which rows?
bByPhyla <- dat$phylum%in% marineNoPlant

#B. classes
needClass <- c(
  "Artiopoda",
  "Branchiopoda",
  "Cephalocarida",
  "Copepoda",
  "Malacostraca",
  "Maxillopoda",
  "Megacheira",
  "Merostomoidea",
  "Ostracoda",
  "Paratrilobita",
  "Pycnogonida",
  "Remipedia",
  "Thylacocephala",
  "Trilobita",
  "Xiphosura"
)

# which rows?
bNeedClass <- dat$class %in% needClass

# subset the data
dat <- dat[bByPhyla | bNeedClass , ]

# remove potential homonymy problems
dat$clgen <- paste(dat$class, dat$genus)

# number of remaining rows
nrow(dat)

# Check to make sure only marine invertebrates are left
unique(dat$phylum)

# remove non-marine environments
omitEnv <- c(
  "\"floodplain\"",
  "alluvial fan",
  "cave",
  "\"channel\"",
  "channel lag" ,
  "coarse channel fill",
  "crater lake",
  "crevasse splay",
  "dry floodplain",
  "delta plain",
  "dune",
  "eolian indet.",
  "fine channel fill",
  "fissure fill",
  "fluvial indet.",
  "fluvial-lacustrine indet.",
  "fluvial-deltaic indet.",
  "glacial",
  "interdune",
  "karst indet.",
  "lacustrine - large",
  "lacustrine - small",
  "lacustrine delta front",
  "lacustrine delta plain",
  "lacustrine deltaic indet.",
  "lacustrine indet.",
  "lacustrine interdistributary bay",
  "lacustrine prodelta",
  "levee",
  "loess",
  "mire/swamp",
  "pond",
  "sinkhole",
  "spring",
  "tar",
  "terrestrial indet.",
  "wet floodplain"
)

# actual omission
dat<-dat[!dat$environment%in%omitEnv, ]

# Filter genera with only one entry
gencount <- as.data.frame(table(dat$genus))

multgen <- gencount$Var1[gencount$Freq > 1]
multgen <- as.character(multgen)

dat <- dat[dat$genus %in% multgen, ]

# Duration per genus

duration <- c()

for (i in unique(dat$genus)) {
  # use upper max age and lower min age to be generous
  maxage <- max(dat$max_ma[dat$genus == i])
  minage <- min(dat$min_ma[dat$genus == i])
  duration[i] <- maxage - minage
}

head(duration)

# Histogram of durations
hist(duration)

# Find the distance between oldest and youngest occurrence

distance <- c()

for (i in unique(dat$genus)) {
  oldlat <- dat$paleolat[dat$genus == i & (dat$max_ma[dat$genus == i] == max(dat$max_ma[dat$genus == i]))]
  oldlong <- dat$paleolng[dat$genus == i & (dat$max_ma[dat$genus == i] == max(dat$max_ma[dat$genus == i]))]
  younglat <- dat$paleolat[dat$genus == i & (dat$min_ma[dat$genus == i] == min(dat$min_ma[dat$genus == i]))]
  younglong <- dat$paleolng[dat$genus == i & (dat$min_ma[dat$genus == i] == max(dat$min_ma[dat$genus == i]))]
  # Find haversine distance
  # Sometimes there are multiple occurrences of same age so only include first entry
  # So that the distance can actually be calculated and not return an error
  point1 <- c(oldlong[1], oldlat[1])
  point2 <- c(younglong[1], younglat[1])
  
  distance[i] <- distHaversine(point1, point2)
  
}

head(distance)

# Histogram of the distances
hist(distance, main = "Genus Distance from Origin at Extinction",
     xlab = "Distance(m)")

results <- cbind(duration, distance)

# log results
logresults <- as.data.frame(log10(results))

# Clean data to remove
# Kinda problematic because it is possible that there are enough 0/infinities to affect results
# The dataset is so big this hopefully won't happen
logresults <- logresults[!is.na(logresults$distance),]
logresults <- logresults[logresults$distance != -Inf, ]

hist(logresults$duration)
hist(logresults$distance)

# Plot Distance and Duration
plot(duration, distance,
     xlab = "Genus Lifespan (Ma)", ylab = "Distance from Location of Genus Origin at Extinction (m)",
     main = "Relationship between Genus Lifespan and Dispersion",
     col = "#1A611E50")

plot(logresults$duration, logresults$distance,
     xlab = expression("log(Genus Lifespan (Ma))"), ylab = "log(Distance from Location of Genus Origin at Extinction (m))",
     main = "Relationship between Genus Lifespan and Dispersion (log)",
     col = "#7C7CD950")

# Test for correlation between duration and distance

cor.test(duration, distance)

cor.test(logresults$duration, logresults$distance)

# Put all the results into a pdf
# Make sure to change path if editing for personal use
pdf("C:/users/jonny/Documents/R projects/Research Project Implementation/Results/graphs.pdf")

hist(duration, main = "Distribution of Genus Durations",
     xlab = "Genus Duration (Ma)")

hist(distance, main = "Genus Distance from Origin at Extinction",
     xlab = "Distance(m)")

plot(duration, distance,
     xlab = "Genus Lifespan (Ma)", ylab = "Distance from Location of Genus Origin at Extinction (m)",
     main = "Relationship between Genus Lifespan and Dispersion",
     col = "#1A611E50")

plot(logresults$duration, logresults$distance,
     xlab = expression("log(Genus Lifespan (Ma))"), ylab = "log(Distance from Location of Genus Origin at Extinction (m))",
     main = "Relationship between Genus Lifespan and Dispersion (log)",
     col = "#7C7CD950")

dev.off()
