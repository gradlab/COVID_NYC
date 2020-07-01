# =============================================================================
# Import packages and key functions
# =============================================================================

# Load key packages
library(tidyverse)
library(lubridate)
library(sf)
library(tidycensus,quietly = T)
library(gridExtra)
library(scales)
library(dtplyr)
library(rmapshaper)
library(broom)
library(purrr)

source('serofunctions2.R')

makesymptomatic <- function(x){
	mutate(x, symptomatic=case_when(primarysx=="None"~0, TRUE~1))
	}

cleanboroughnames <- function(x){
	x %>% mutate(borough=case_when(
		borough=="NorthQueens"~"North Queens",
		borough=="SouthQueens"~"South Queens",
		borough=="StatenIsland"~"Staten Island",
		TRUE~borough
		)) %>%
	mutate(borough=factor(borough, levels=c("Bronx","Brooklyn","Manhattan","North Queens","South Queens","Overall")))
}

restrictdates <- function(x){
	x %>%
	filter(date>=ymd("2020-03-16")) %>%
	filter(date<ymd("2020-04-20"))
}

# Load in personal Census API key (new users will need to apply for their own and define the censusapi variable here) 
source("data/censusapi.R") 

# =============================================================================
# Import L&D COVID data
# =============================================================================

dat_agg <- read.csv("data/nyc_hospital_data.csv", stringsAsFactors=FALSE) %>% as_tibble()
dat_agg$firstdate <- as_date(dat_agg$firstdate) 

# =============================================================================
# Import movement data
# =============================================================================

mvmtdat <- data.frame(borough=c("SouthQueens","Bronx","Brooklyn","NorthQueens","Manhattan"), perc_change=c(0.414,0.495,0.524,0.574,0.687), stringsAsFactors=FALSE) 

# =============================================================================
# Import mapping data
# =============================================================================

census_api_key(censusapi)
options(tigris_use_cache = TRUE)

# Import ZIP boundaries:
ny_zips_raw <- get_acs(geography = "zcta", geometry = TRUE, variables = "B01003_001") %>% st_transform(crs=3857)
# Import tracts to dissolve into state boundary
ny_tracts <- get_acs(state = "NY", geography = "tract", geometry = TRUE, variables = "B01003_001") %>% st_transform(crs=3857)
# Dissolve tracts into state boundary
ny_state <- ms_dissolve(ny_tracts)
# Keep only the ZIPs in New York state:
ny_zips <- st_intersection(ny_state, ny_zips_raw)

# Import mapping between boroughs and zips:
boroughzip <- read.csv("data/borough_zip.csv", colClasses=c("character","character"))
# Import mapping between fine boroughs (including N/S Queens) and 3-digit zips:
boroughzip3 <- read.csv("data/borough_zip3.csv", colClasses=c("character","character"))
# Import mapping between fine boroughs (including N/S Queens) and 3-digit zips:
boroughzip3fine <- read.csv("data/borough_zip3_fine.csv", colClasses=c("character","character"))

# Generate geometries for the fine boroughs: 
ny_boroughlist <- ny_zips %>% 
	right_join(boroughzip, by=c("GEOID"="zip")) %>%
	filter(st_is(geometry, c("POLYGON", "MULTIPOLYGON"))) %>%
	mutate(zip3=substr(GEOID,1,3)) %>%
	select(-borough) %>%
	inner_join(boroughzip3, by="zip3") %>%
	split(.$borough) %>%
	map(~ ms_dissolve(.)) 
ny_borough <- ny_boroughlist %>%
	reduce(rbind) %>%
	mutate(borough=names(ny_boroughlist)) %>%
	select(-rmapshaperid)

# Generate geometries for the fine boroughs: 
ny_boroughlist_fine <- ny_zips %>% 
	right_join(boroughzip, by=c("GEOID"="zip")) %>%
	filter(st_is(geometry, c("POLYGON", "MULTIPOLYGON"))) %>%
	mutate(zip3=substr(GEOID,1,3)) %>%
	select(-borough) %>%
	inner_join(boroughzip3fine, by="zip3") %>%
	split(.$borough) %>%
	map(~ ms_dissolve(.)) 
ny_borough_fine <- ny_boroughlist_fine %>%
	reduce(rbind) %>%
	mutate(borough=names(ny_boroughlist_fine)) %>%
	select(-rmapshaperid)

# Save the bounding box around these ZIPs:
boroughbox <- st_bbox(ny_borough_fine)

# Store the hospital coordinates: 
hosplocs <- data.frame(site=c("NYP-WCM","NYP-Queens", "NYP-LMH", "NYP-CUIMC","MSH","MSW"), 
	geometry=c(
		st_sfc(st_point(c(-73.951953, 40.768786)), crs=4326),
		st_sfc(st_point(c(-73.825319, 40.747239)), crs=4326),
		st_sfc(st_point(c(-74.004769, 40.709901)), crs=4326),
		st_sfc(st_point(c(-73.941967, 40.841752)), crs=4326),
		st_sfc(st_point(c(-73.952508, 40.789940)), crs=4326),
		st_sfc(st_point(c(-73.986497, 40.769784)), crs=4326)
		)) %>%
	st_sf() %>%
	st_transform(crs=3857)

hosplocsdf <- cbind(st_drop_geometry(hosplocs), 
	do.call(rbind, st_geometry(hosplocs)) %>% 
    as_tibble() %>% 
    setNames(c("lon","lat")))

labelhosp <- list(
	geom_sf(data=hosplocs, aes(geometry=geometry, fill=NULL)),
	geom_label(data=hosplocsdf, aes(x=lon, y=lat, label=site, geometry=NULL), vjust=c(0,0,0,0,0,0), hjust=c(0,0,1,1,1,1), nudge_x=c(1000, 1000, -1000, -1000,-1000,-1000), alpha=0.6, size = 3,
                   label.padding = unit(0.1, "lines"), fill="white")
	)

# =============================================================================
# Basic statistics
# =============================================================================

# Number of women
dat_agg %>% summarise(sum(N))

# Covid result
dat_agg %>% group_by(covidresult) %>% summarise(sum(N))

# Borough 
dat_agg %>% group_by(borough) %>% summarise(sum(N))

# =============================================================================
# Fig.1: posterior prevalence
# =============================================================================

se <- 0.9
sp <- 1

# Overall: 
postpopprev_borough_fine <- dat_agg %>% 
	group_by(borough, covidresult) %>%
	summarise(N=sum(N)) %>%
	ungroup() %>%
	mutate(covidresult=case_when(covidresult==1~"pos", TRUE~"neg")) %>%
	pivot_wider(names_from=covidresult, values_from=N) %>%
	rbind(data.frame(borough="Overall",pos=sum(.$pos),neg=sum(.$neg))) %>%
	split(.$borough) %>% 
	map(~ data.frame(post=sample_post_r_log(.$pos,.$neg,se,sp,10000))) %>%
	bind_rows(.id="borough")

postpopprev_borough_fine_summ <- postpopprev_borough_fine %>% 
	group_by(borough) %>%
	summarise(postmean=mean(post), postlwr=quantile(post, 0.025), postupr=quantile(post, 0.975))

figpostpopprev_borough_fine <- postpopprev_borough_fine %>% 
	cleanboroughnames %>%
	rename("Borough"=borough) %>%
	mutate(post=100*post) %>%
	ggplot(aes(x=post, col=Borough)) + 
		# geom_density() + 
		stat_density(geom="line",position="identity", adjust=2) + 
		xlim(c(0,75)) + 
		labs(x="Estimated prevalence (%)", y="Probability density") +
		scale_colour_manual(values=c("Bronx"="green","Brooklyn"="blue","Manhattan"="magenta","South Queens"="red","North Queens"="purple","Staten Island"="brown","Overall"="black")) + 
		theme_minimal() + 
		theme(text=element_text(size=14))

ggsave(figpostpopprev_borough_fine,file="figures/postpopprev_borough_fine.pdf", width=7, height=4)

# Weekly: 

postpopprev_borough_fine_wk <- dat_agg %>% 
	group_by(borough, firstdate) %>% 
	mutate(covidresult=case_when(covidresult==1~"pos", TRUE~"neg")) %>%
	pivot_wider(names_from=covidresult, values_from=N) %>%
	replace_na(list(neg=0, pos=0)) %>%
	group_by(firstdate) %>%
	bind_rows(summarise(., borough="Overall",pos=sum(pos),neg=sum(neg))) %>%
	split(.$borough) %>%
	map(~ split(., .$firstdate)) %>%
	map_depth(2, ~ data.frame(post=sample_post_r_log(.$pos,.$neg,se,sp,10000))) %>%
	map_depth(2, ~ data.frame(
		mean=mean(.$post), 
		lwr=quantile(.$post, 0.025),
		upr=quantile(.$post, 0.975))) %>%
	map(~ bind_rows(., .id="firstdate")) %>%
	bind_rows(.id="borough") %>%
	mutate(firstdate = ymd(firstdate))

figpostpopprev_borough_fine_wk <- postpopprev_borough_fine_wk %>% 
	cleanboroughnames %>%
	ggplot(aes(x=firstdate, y=100*mean, col=borough)) + 
		geom_errorbar(aes(ymin=100*lwr, ymax=100*upr), width=1.5, alpha=0.5) + 
		geom_point() + 
		geom_line() + 
		ylim(c(0,100)) + 
		scale_x_date(breaks="1 week", minor_breaks="1 week", date_labels = "%b %d") + 
		scale_colour_manual(values=c("Bronx"="green","Brooklyn"="blue","Manhattan"="magenta","South Queens"="red","North Queens"="purple","Staten Island"="brown","Overall"="black")) + 
		theme_minimal() + 
		theme(text=element_text(size=16), axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5),legend.position="none") + 
		labs(x="Week", y="Estimated prevalence (%)") + 
		facet_wrap(vars(borough), ncol=2)

ggsave(figpostpopprev_borough_fine_wk,file="figures/postpopprev_borough_fine_wk.pdf", width=7, height=8)


# =============================================================================
# Fig2: maps 
# =============================================================================


figprevmap <- postpopprev_borough_fine_summ %>%
	mutate(postmean=100*postmean) %>%
	left_join(ny_borough_fine) %>%
	rename("Estimated\nprevalence\n(%)"=postmean) %>%
	ggplot(aes(geometry=geometry, fill=`Estimated\nprevalence\n(%)`)) + 
		geom_sf(data=ny_state, fill="lightgray", color="gray") + 
		geom_sf() + 
		labelhosp + 
		scale_fill_distiller(palette="Blues", direction=1) + 
		coord_sf(xlim=c(boroughbox$xmin, boroughbox$xmax), ylim=c(boroughbox$ymin, boroughbox$ymax)) + 
		theme_void()

ggsave(figprevmap,file="figures/prevmap.pdf")


# The updated cumulative flux: 
figcummvmt <- mvmtdat %>%
	mutate(perc_change=100*perc_change) %>%
	select(borough, perc_change) %>%
	left_join(ny_borough_fine) %>%
	rename("Percent decline in\ncommuting movements"=perc_change) %>%
	ggplot(aes(geometry=geometry, fill=`Percent decline in\ncommuting movements`)) + 
		geom_sf(data=ny_state, fill="lightgray", color="gray") + 
		geom_sf() + 
		labelhosp + 
		scale_fill_distiller(palette="Blues", direction=-1) + 
		coord_sf(xlim=c(boroughbox$xmin, boroughbox$xmax), ylim=c(boroughbox$ymin, boroughbox$ymax)) + 
		theme_void() 
ggsave(figcummvmt,file="figures/cummvmt.pdf")
	

# =============================================================================
# TableS2: Weekly tests and positives by week and borough
# =============================================================================

fullfirstdates <- data.frame(date=seq(ymd("2020-03-16"), ymd("2020-04-27"), by="1 week"))

dat_agg %>%
	group_by(borough, firstdate) %>% 
	summarise(pos=sum(N*covidresult), neg=sum(-N*(covidresult-1)), ntests=sum(N), ppos=pos/ntests) %>%
	group_by(firstdate) %>% 
	bind_rows(summarise(., borough="Overall",pos=sum(pos, na.rm=TRUE),neg=sum(neg, na.rm=TRUE), ntests=sum(ntests,na.rm=TRUE), ppos=pos/ntests)) %>% 
	ungroup() %>% 
	split(.$borough) %>%
	map(~ select(.,-borough)) %>%
	map(~ right_join(., fullfirstdates, by=c("firstdate"="date"))) %>%
	map(~ mutate(., ppos=round(ppos*1000)/10)) %>%
	map(~ replace_na(., list(pos="-", neg="-", ntests="-", ppos="-"))) %>%
	map(~ select(., firstdate, ntests, pos, ppos)) %>%
	bind_rows(.id="borough") %>%
	mutate(firstdate=stamp("March 1")(firstdate)) %>%
	cleanboroughnames %>%
	arrange(borough) %>%
	print(n=50)


# =============================================================================
# TableS3: posterior prev by borough w/ varying sensitivity/specificity
# =============================================================================

sp <- 1

postpopprev_borough_fine_multise <- dat_agg %>% 
	group_by(borough) %>% 
	summarise(pos=sum(N*covidresult), neg=sum(-N*(covidresult-1))) %>%
	rbind(data.frame(borough="Overall", pos=sum(.$pos), neg=sum(.$neg))) %>%
	split(.$borough) %>% 
	map(~ data.frame(
		post70=sample_post_r_log(.$pos,.$neg,0.7,sp,10000),
		post80=sample_post_r_log(.$pos,.$neg,0.8,sp,10000),
		post90=sample_post_r_log(.$pos,.$neg,0.9,sp,10000)
		)) %>%
	bind_rows(.id="borough") %>%
	pivot_longer(c("post70","post80","post90"), names_to="sensitivity", values_to="post") %>%
	mutate(sensitivity = substr(sensitivity, 5,6))

postpopprev_borough_fine_multise %>%	
	cleanboroughnames %>%
	split(.$sensitivity) %>%
	map(~ group_by(., borough)) %>%
	map(~ summarise(., postmean=100*mean(post), postlwr=100*quantile(post, 0.025), postupr=100*quantile(post, 0.975))) %>%
	bind_rows(.id="sensitivity")

# =============================================================================
# FigS1: Regression
# =============================================================================

figregression <- postpopprev_borough_fine_summ %>%
	filter(borough!="Overall") %>%
	left_join(mvmtdat, by="borough") %>%
	mutate(postmean=100*postmean) %>%
	mutate(postlwr=100*postlwr) %>%
	mutate(postupr=100*postupr) %>%
	mutate(perc_change=100*perc_change) %>%
	ggplot(aes(x=perc_change, y=postmean)) + 
		geom_point() + 
		geom_line(stat="smooth", method="lm") + 
		geom_errorbar(aes(ymin=postlwr, ymax=postupr), width=1.5, alpha=0.5) + 
		labs(x="Percent decline in commuting movements", y="Estimated prevalence (%)") + 
		theme_minimal() + 
		theme(text=element_text(size=18))
ggsave(figregression,file="figures/regression.pdf")

temp <- postpopprev_borough_fine_summ %>%
	filter(borough!="Overall") %>%
	left_join(mvmtdat, by="borough") %>%
	mutate(postmean=100*postmean) %>%
	mutate(postlwr=100*postlwr) %>%
	mutate(postupr=100*postupr) %>%
	mutate(perc_change=100*perc_change)

# Generate regression coefficient with uncertainty: 

corrcoeffs <- postpopprev_borough_fine %>%
	filter(borough!="Overall") %>%
	group_by(borough) %>%
	sample_n(10000) %>%
	left_join(mvmtdat, by="borough") %>%
	group_by(borough) %>%
	mutate(row=1:n()) %>%
	split(.$row) %>%
	map(~ data.frame(R=cor(.$post, .$perc_change))) %>%
	bind_rows(.id="row") %>%
	as_tibble()

figRdensity <- corrcoeffs %>% 
	ggplot(aes(x=R)) + 
		geom_density(adjust=3) + 
		labs(x="Regression coefficient (R)", y="Probability density") + 
		theme_minimal() + 
		theme(axis.ticks.y=element_blank(), axis.text.y=element_blank())

corrcoeffs %>%
	summarise(Rmean=mean(R), Rlwr=quantile(R, 0.025), Rupr=quantile(R, 0.975))

