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
results <- as.data.frame(results)
results$genus <- unique(dat$genus)

# log results
logresults <- as.data.frame(log10(results))

# Clean data to remove
# Kinda problematic because it is possible that there are enough 0/infinities to affect results
# The dataset is so big this hopefully won't happen
logresults <- logresults[!is.na(logresults$distance),]
logresults <- logresults[logresults$distance != -Inf, ]

hist(logresults$duration)
hist(logresults$distance)

# Statistics of Results
resultstat <- lm(duration ~ distance, results)
summary(resultstat)

logresultstat <- lm(duration ~ distance, logresults)
summary(logresultstat)

# Plot Distance and Duration
plot(distance, duration,
     ylab = "Genus Lifespan (Ma)", xlab = "Distance from Location of Genus Origin at Extinction (m)",
     main = "Relationship between Genus Lifespan and Dispersion",
     col = "#1A611E50")

plot(y = logresults$duration, x = logresults$distance,
     ylab = expression("log(Genus Lifespan (Ma))"), xlab = "log(Distance from Location of Genus Origin at Extinction (m))",
     main = "Relationship between Genus Lifespan and Dispersion (log)",
     col = "#7C7CD950")

# Test for correlation between duration and distance

cor.test(duration, distance)

cor.test(logresults$duration, logresults$distance)




####################
# Compare different groups

# Comparing tropical to non-tropical taxa
troplat <- 23.4 # establishing tropics at 23.4 degrees latitude

tropicdat <- abs(dat$paleolat) <= troplat

nontropicgenera <- unique(dat$genus[!tropicdat]) # Species which do not rely on the tropics but may enter it
tropicgenera <- unique(dat$genus[!dat$genus %in% nontropicgenera]) # Species which rely on tropics

# Compare duration/distances for tropical/nontropical genera
tropicresults <- results[results$genus %in% tropicgenera, ]

nontropicresults <- results[results$genus %in% nontropicgenera, ]

logtropicresults <- tropicresults
logtropicresults$duration <- log10(logtropicresults$duration)
logtropicresults$distance <- log10(logtropicresults$distance)
logtropicresults <- logtropicresults[!is.na(logtropicresults),]
logtropicresults <- logtropicresults[logtropicresults$distance != -Inf,]

lognontropicresults <- nontropicresults
lognontropicresults$duration <- log10(lognontropicresults$duration)
lognontropicresults$distance <- log10(lognontropicresults$distance)
lognontropicresults <- lognontropicresults[!is.na(lognontropicresults),]
lognontropicresults <- lognontropicresults[lognontropicresults$distance != -Inf,]

# Histograms of duration and distance for tropic and nontropic genera
hist(tropicresults$duration)
hist(tropicresults$distance)

hist(logtropicresults$duration)
hist(logtropicresults$distance)

hist(nontropicresults$duration)
hist(nontropicresults$distance)

hist(lognontropicresults$duration)
hist(lognontropicresults$distance)

# Statistics of tropic results
tropicstat <- lm(duration ~ distance, logtropicresults)
summary(tropicstat)

nontropicstat <- lm(duration ~ distance, lognontropicresults)
summary(nontropicstat)

# Scatterplots for results
plot(tropicresults$duration, tropicresults$distance)
plot(y = logtropicresults$duration, x = logtropicresults$distance,
     main = "Dispersion of Tropical Genera",
     ylab = "log(Genus Duration(Ma))",
     xlab = "Log(Distance from Origin at Extinction(m))",
     col = "#34eb3480", pch = 16)
abline(tropicstat)

plot(nontropicresults$duration, nontropicresults$distance)
plot(y = lognontropicresults$duration, x = lognontropicresults$distance,
     main = "Dispersion of Non-Tropical Genera",
     ylab = "log(Genus Duration(Ma))",
     xlab = "Log(Distance from Origin at Extinction(m))",
     col = "#38e8e850", pch = 16)
abline(nontropicstat)

################################
# Compare results for taxa which are regarded as disperers vs nondispersers
# Using lecithotrophic vs planktotrophic brachiopods
plankorder <- c("Lingulida", "Discinida")

plankbrach <- dat$genus[dat$order %in% plankorder]

lecithorder <- c("Craniida", "Terebratulida", "Rhynchonellida")

lecithbrack <- dat$genus[dat$order %in% lecithorder]

# Get logged stats for both  groups
plankresults <- results[results$genus %in% plankbrach, ]
plankresults$duration <- log10(plankresults$duration)
plankresults$distance <- log10(plankresults$distance)
plankresults <- plankresults[plankresults$distance != -Inf, ]
plankstats <- lm(duration ~ distance, plankresults)
summary(plankstats)

lecitresults <- results[results$genus %in% lecithbrack, ]
lecitresults$duration <- log10(lecitresults$duration)
lecitresults$distance <- log10(lecitresults$distance)
lecitresults <- lecitresults[lecitresults$distance != -Inf, ]
lecitstats <- lm(duration ~ distance, lecitresults)
summary(lecitstats)

# Plot results
plot(plankresults$distance, plankresults$duration,
     main = "Duration of Planktotrophic Brachiopods",
     xlab = "log(Distance from Origin at Extinction(m))",
     ylab = "log(Duration of Genus(Ma))",
     col = "#0004fc", pch = 16)
abline(plankstats)

plot(lecitresults$distance, lecitresults$duration,
     main = "Duration of Lecithotrophic Brachiopods",
     xlab = "log(Distance from Origin at Extinction(m))",
     ylab = "log(Duration of Genus(Ma))",
     col = "#fc5000", pch = 16)
abline(lecitstats)

###########################################
# Put all the results into a pdf
# Make sure to change path if editing for personal use
pdf("C:/users/jonny/Documents/R projects/Research Project Implementation/Results/graphs.pdf")

hist(duration, main = "Distribution of Genus Durations",
     xlab = "Genus Duration (Ma)")

hist(distance, main = "Genus Distance from Origin at Extinction",
     xlab = "Distance(m)")

plot(distance, duration,
     ylab = "Genus Lifespan (Ma)", xlab = "Distance from Location of Genus Origin at Extinction (m)",
     main = "Relationship between Genus Lifespan and Dispersion",
     col = "#1A611E50")
abline(resultstat)
print(summary(resultstat))

plot(y = logresults$duration, x = logresults$distance,
     ylab = expression("log(Genus Lifespan (Ma))"), xlab = "log(Distance from Location of Genus Origin at Extinction (m))",
     main = "Relationship between Genus Lifespan and Dispersion (log)",
     col = "#7C7CD950")
abline(logresultstat)
summary(logresultstat)

plot(y = logtropicresults$duration, x = logtropicresults$distance,
     main = "Dispersion of Tropical Genera",
     ylab = "log(Genus Duration(Ma))",
     xlab = "Log(Distance from Origin at Extinction(m))",
     col = "#34eb3480", pch = 16)
abline(tropicstat)
print(summary(tropicstat))

plot(y = lognontropicresults$duration, x = lognontropicresults$distance,
     main = "Dispersion of Non-Tropical Genera",
     ylab = "log(Genus Duration(Ma))",
     xlab = "Log(Distance from Origin at Extinction(m))",
     col = "#38e8e850", pch = 16)
abline(nontropicstat)
print(summary(nontropicstat))

plot(plankresults$distance, plankresults$duration,
     main = "Duration of Planktotrophic Brachiopods",
     xlab = "log(Distance from Origin at Extinction(m))",
     ylab = "log(Duration of Genus(Ma))",
     col = "#0004fc", pch = 16)
abline(plankstats)
print(summary(plankstats))

plot(lecitresults$distance, lecitresults$duration,
     main = "Duration of Lecithotrophic Brachiopods",
     xlab = "log(Distance from Origin at Extinction(m))",
     ylab = "log(Duration of Genus(Ma))",
     col = "#fc5000", pch = 16)
abline(lecitstats)
print(summary(lecitstats))

dev.off()
