################ Analysis for Genomic Sentinile site selection #######

source(here::here("scripts", "0_functions.R"))
pacman::p_load(
  rio,           # to import data
  here,          # to locate files
  tidyverse,     # to clean, handle, and plot the data (includes ggplot2 package)
  sf,            # to manage spatial data using a Simple Feature format
  tmap,          # to produce simple maps, works for both interactive and static maps
  janitor,       # to clean column names
  #OpenStreetMap, # to add OSM basemap in ggplot map
  spdep ,         # spatial statistics,
  here,
  geojsonio
)

dat <- load_clean_data()
detach("package:import", unload = TRUE)
#Extracting the data to work on in to one data frame
# The disease include Malaria, Anthrax, Meningitis, afp_polio,rabies(exposed),
#Cholera, Measles and Dengue Fever

msls <-  import(here("raw data/Linelist/Measles_linelsit_aggregate.csv"))
chol <- import(here("raw data/Linelist/Cholera22_24_linelsit_aggregate.csv"))
den <- import(here("raw data/Linelist/Dengue_linelsit_aggregate_since2023.csv"))

# Start with the malaria data and add the rest
dat_gen <- dat %>% select(region,zone,  woreda, ADM3_PCODE, year, epi_week,
                          totl_malr, mnng_totl,afp_poli,anth,rabs
                          
) %>% 
  filter(year >= 2020 & year < 2025 ) %>%  group_by(region, zone, woreda, ADM3_PCODE, year) %>% 
  summarise(mal = sum(totl_malr, na.rm = T),
            mng = sum(mnng_totl, na.rm = T),
            afp = sum(afp_poli, na.rm = T),
            anth = sum(anth, na.rm = T),
            rab = sum(rabs, na.rm = T), .groups = "drop")

# pop <- import(here("Pop/Population by GIS woreda.xlsx"))

# dat_2024 <- dat_2024 %>% left_join(pop, by =c("woreda" = "Gnames")) %>% janitor::clean_names() 
# sum(is.na(dat_2024$population))
# dat_2024 <- dat_2024 %>% mutate(mal_api = round((mal/population) * 1000, 1))

#Filter out date from 2020 to 2024 data 
msls <- msls %>% filter(year >= 2020 & year < 2025) %>% 
  group_by(region, zone, woreda, ADM3_PCODE, year) %>% 
  summarise(msls = sum(msls, na.rm = T),
            .groups = "drop") 
chol <- chol %>% filter(year_dr >= 2020 & year_dr < 2025) %>%  #but only year 22 to 24 available
  group_by(region, zone, woreda, ADM3_PCODE, year_dr) %>% 
  summarise(chlr = sum(chlr, na.rm = T),
            .groups = "drop") %>% 
  dplyr::rename(year = year_dr)

den <- den %>%  filter(year >= 2020 & year < 2025) %>% # but only avaialble from 23 to 24
  group_by(region, zone, woreda, ADM3_PCODE, year) %>% 
  summarise(deng_fevr = sum(deng_fevr, na.rm = T),
            .groups = "drop")

# To quickly see the sum of dengue cases in the main data frame
sum(dat[dat$year == 2024, "deng_fevr"],na.rm = T)
sum(dat[dat$year == 2024, "deng_fevr_dths"],na.rm = T)
summary(dat[dat$year == 2024, "deng_fevr"],na.rm = T)

# Now moving on with merging the data sets


dat_gen <- dat_gen  %>% 
  full_join(msls, by = c("region", "zone", "woreda", "ADM3_PCODE", "year"))  
dat_gen <- dat_gen %>% full_join(chol, by = c("region", "zone", "woreda", "ADM3_PCODE", "year"))
dat_gen <- dat_gen %>% full_join(den, by = c("region", "zone", "woreda", "ADM3_PCODE", "year"))

glimpse(dat_gen)

# Proceeding with the analysis - identifying top burdne woredas 
dat_agg <- dat_gen %>%
  group_by(region, zone, woreda, ADM3_PCODE) %>%
  summarise(
    mal = sum(mal, na.rm = TRUE),
    mng = sum(mng, na.rm = TRUE),
    afp = sum(afp, na.rm = TRUE),
    anth = sum(anth, na.rm = TRUE),
    rab = sum(rab, na.rm = TRUE),
    msls = sum(msls, na.rm = TRUE),
    chlr = sum(chlr, na.rm = TRUE),
    deng_fevr = sum(deng_fevr, na.rm = TRUE),
    .groups = "drop"
  )

#Identifying Top woredas for the multiple disesaes
disease_cols <- c("mal", "mng", "afp", "anth", "rab", "msls", "chlr", "deng_fevr")
top_woredas <- list()
for (disease in disease_cols) {
  top_woredas[[disease]] <- dat_agg %>%
    select(region, zone, woreda, ADM3_PCODE, !!sym(disease)) %>%
    filter(!!sym(disease) > 0) %>%
    group_by(region) %>%
    arrange(desc(!!sym(disease))) %>%
    mutate(rank = row_number()) %>%
    filter(rank <= if_else(region %in% c("Oromia", "Amhara"), 20L, 10L)) %>%
    mutate(disease = disease) %>%
    select(region, zone, woreda, ADM3_PCODE, disease, !!sym(disease), rank)
}

top_woredas_df <- bind_rows(top_woredas)
# To export the data by each woreda
export(top_woredas_df %>% ungroup() %>% relocate(disease, .before = region), here("output/Top burden woredas by region.csv"))


dat_multi <- top_woredas_df %>%
  group_by(region, zone, woreda, ADM3_PCODE) %>%
  summarise(
    num_diseases = n_distinct(disease),
    diseases_present = paste(disease, collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(num_diseases >= 2) %>%
  #To join with the aggregate data and sum all the disease 
  left_join(
    dat_agg %>% 
      mutate( total_cases = rowSums(select(.,all_of(disease_cols)),na.rm = T)) %>% 
      select(region,zone, woreda, ADM3_PCODE, total_cases),
    by = c("region", "zone", "woreda", "ADM3_PCODE")
  ) %>% 
  arrange(region, desc(num_diseases), desc(total_cases))

#Sites with 2 or more disease in each region prioritized based on the 
#number of diseases and summation of all disease cases
export( dat_multi, here("output/Top burden woredas by region with disease category and total case.csv"))

#Working on a sumamry map
eth_adm3_raw <- read_sf(here("shapefile/eth_admbnda_adm3_csa_bofedb_2021.shp"))
eth_adm3 <- eth_adm3_raw %>% janitor::clean_names()

dat_multi_sf <- eth_adm3 %>%
  left_join(dat_multi, by = c("adm3_pcode" = "ADM3_PCODE")) %>%
  filter(!is.na(num_diseases)) %>% # Keep only woredas with data
  st_as_sf()

tmap_mode("view")

sentinel_map <- tm_shape(dat_multi_sf) +
  tm_polygons(
    "num_diseases",
    title = "Number of Diseases",
    palette = "Blues",
    popup.vars = c("region", "zone", "woreda", "num_diseases", "diseases_present", "total_cases"),
    id = "woreda"
  ) +
  tm_borders() +
  tm_layout(
    title = "Sentinel Sites with Multiple Diseases (2020-2024)",
    title.position = c("center", "top")
  )
sentinel_map

eth_adm1_raw <- read_sf(here("shapefile/eth_admbnda_adm1_csa_bofedb_2021.shp"))
eth_adm1 <- eth_adm1_raw %>% janitor::clean_names()

# Pre-compute centroids for Admin 1 labels
eth_adm1_centroids <- eth_adm1 %>%
  mutate(
    centroid = st_centroid(geometry),
    lng = st_coordinates(centroid)[,1],
    lat = st_coordinates(centroid)[,2]
  ) %>%
  select(adm1_en, lng, lat)

sentinel_map2 <- tm_shape(dat_multi_sf) +
  tm_polygons(
    "num_diseases",
    title = "Number of Diseases",
    palette = "Blues",
    popup.vars = c("region", "zone", "woreda", "num_diseases", "diseases_present", "total_cases"),
    id = "woreda",
    alpha = 0.7
  ) +
  tm_shape(eth_adm1) +
  tm_borders(
    col = "black",
    lwd = 2,
    alpha = 0.8
  ) +
  tm_text(
    "adm1_en",
    size = 0.8,
    col = "black",
    fontface = "bold"
  ) +
  tm_layout(
    title = "Sentinel Sites with Multiple Diseases (2020-2024)",
    title.position = c("center", "top"),
    # legend.position = c("left", "bottom")
  )
sentinel_map2

tmap_save(sentinel_map2, filename = here("output/sentinel_map_multi 2.html"))

#Saving processed data for qmd
saveRDS(dat_multi, file = here("output/dat_multi.rds"))
saveRDS(dat_multi_sf, file = here("output/dat_multi_sf.rds"))
saveRDS(eth_adm1, file = here("output/eth_adm1.rds"))


sentinel_map2





# Trying to add further functionality like selecting the number of disease ######
# Load required packages
pacman::p_load(
  tidyverse,      # Data manipulation and visualization
  sf,             # Spatial data handling
  leaflet,        # Interactive mapping with filtering
  htmlwidgets,    # For custom JS in Leaflet
  htmltools,      # For HTML elements in Leaflet
  DT,             # Interactive tables in QMD
  here,            # File path management
  geojsonio
)
install.packages("geojsonio")
#Step 5: Create an interactive Leaflet map with dropdown filter
# Define color palette for num_diseases (mimicking tmap's "Blues")
pal <- colorNumeric(
  palette = c("#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C"),
  domain = dat_multi_sf$num_diseases
)

# Unique number of diseases for dropdown
disease_counts <- sort(unique(dat_multi_sf$num_diseases))


# Create Leaflet map
leaflet_map <- leaflet(data = dat_multi_sf) %>%
  addTiles() %>%
  addPolygons(
    fillColor = ~pal(num_diseases),
    fillOpacity = 0.7,
    color = "#BDBDC3",
    weight = 1,
    popup = ~paste0(
      "<b>Woreda:</b> ", woreda, "<br>",
      "<b>Region:</b> ", region, "<br>",
      "<b>Zone:</b> ", zone, "<br>",
      "<b>Number of Diseases:</b> ", num_diseases, "<br>",
      "<b>Diseases:</b> ", diseases_present, "<br>",
      "<b>Total Cases:</b> ", total_cases
    ),
    group = "Woredas",
    layerId = ~woreda
  ) %>%
  addPolygons(
    data = eth_adm1,
    fill = FALSE,
    color = "black",
    weight = 2,
    opacity = 0.8,
    group = "Regions"
  ) %>%
  addLabelOnlyMarkers(
    data = eth_adm1_centroids,
    lng = ~lng,
    lat = ~lat,
    label = ~adm1_en,
    labelOptions = labelOptions(
      noHide = TRUE,
      textsize = "12px",
      textOnly = TRUE,
      style = list("font-weight" = "bold")
    ),
    group = "Region Labels"
  ) %>%
  addLegend(
    position = "bottomleft",
    pal = pal,
    values = ~num_diseases,
    title = "Number of Diseases",
    opacity = 0.7
  )

# Add dropdown control using custom HTML and JavaScript
dropdown_html <- tags$div(
  style = "position: absolute; top: 10px; left: 50px; z-index: 1000; background: white; padding: 5px; border-radius: 5px;",
  tags$label("Select Number of Diseases: "),
  HTML(
    sprintf(
      '<select id="diseaseFilter" onchange="filterMap(this.value)">
         <option value="all">All</option>
         %s
       </select>',
      paste(
        sprintf('<option value="%s">%s Diseases</option>', disease_counts, disease_counts),
        collapse = ""
      )
    )
  )
)

# JavaScript for filtering
js_code <- HTML(
  sprintf(
    '
    <script>
    function filterMap(value) {
      var map = window.leafletMap;
      map.eachLayer(function(layer) {
        if (layer.options.group === "Woredas") {
          map.removeLayer(layer);
        }
      });
      var data = %s;
      var filteredData = value === "all" ? data : data.filter(function(d) {
        return d.properties.num_diseases == value;
      });
      L.geoJSON(filteredData, {
        style: function(feature) {
          return {
            fillColor: getColor(feature.properties.num_diseases),
            fillOpacity: 0.7,
            weight: 1,
            color: "#BDBDC3"
          };
        },
        onEachFeature: function(feature, layer) {
          layer.bindPopup(
            "<b>Woreda:</b> " + feature.properties.woreda + "<br>" +
            "<b>Region:</b> " + feature.properties.region + "<br>" +
            "<b>Zone:</b> " + feature.properties.zone + "<br>" +
            "<b>Number of Diseases:</b> " + feature.properties.num_diseases + "<br>" +
            "<b>Diseases:</b> " + feature.properties.diseases_present + "<br>" +
            "<b>Total Cases:</b> " + feature.properties.total_cases
          );
          layer.options.group = "Woredas";
        }
      }).addTo(map);
    }

    function getColor(num_diseases) {
      var colors = %s;
      var breaks = %s;
      for (var i = 0; i < breaks.length; i++) {
        if (num_diseases <= breaks[i]) {
          return colors[i];
        }
      }
      return colors[colors.length - 1];
    }

    document.addEventListener("DOMContentLoaded", function() {
      window.leafletMap = window.leafletMap || document.querySelector(".leaflet")._leaflet_map;
      filterMap("all");
    });
    </script>
    ',
    geojsonio::geojson_json(dat_multi_sf, geometry = "geometry"),
    jsonlite::toJSON(c("#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C")),
    jsonlite::toJSON(sort(unique(dat_multi_sf$num_diseases)))
  )
)

# Add dropdown and JS as dependencies to the widget
leaflet_map <- htmlwidgets::prependContent(
  leaflet_map,
  dropdown_html,
  js_code
)
leaflet_map
# Save the map as an HTML file
saveWidget(leaflet_map, file = here("output/sentinel_map_multi.html"), selfcontained = TRUE)

# Save processed data for QMD
saveRDS(dat_multi, file = here("output/dat_multi.rds"))
saveRDS(dat_multi_sf, file = here("output/dat_multi_sf.rds"))
saveRDS(eth_adm1, file = here("output/eth_adm1.rds"))
saveRDS(eth_adm1_centroids, file = here("output/eth_adm1_centroids.rds"))
Chang




#  Plotting sentinle site maps --------------------------------------------

eth_adm3_raw <- read_sf(here("shapefile/eth_admbnda_adm3_csa_bofedb_2021.shp"))
eth_adm3 <- eth_adm3_raw %>% janitor::clean_names()

eth_adm1 <- read_sf(here("shapefile/eth_admbnda_adm1_csa_bofedb_2021.shp")) %>%  
  janitor::clean_names()

sent_site <- import(here("Senti sites.csv"))

sent_site_sf <- sent_site %>% filter(sentinile == "yes") %>% 
  left_join(eth_adm3, by = c("ADM3_PCODE" = "adm3_pcode")) %>%
  select(region, zone, woreda, ADM3_PCODE, sentinile, diseases, laboratory, geometry)%>% # Keep only woredas with data
  st_as_sf()

tmap_mode("plot")


senti_map <- tm_shape(sent_site_sf) +
  
  # # faint ADM3 boundaries underneath
  # #tm_borders(
  #   lwd       = 0.3,
  #   col       = "grey80",
  #   layer.name= "ADM3 boundaries"
  # ) +
  
  # fill choropleth of cholera cases with fixed breaks
  tm_polygons(
    "sentinile",
    palette    = "green", #magma , inferno 
    fill.scale = tm_scale_categorical(),
    lwd = 0.3,
    #title      = "Cholera case",
    border.col = NA,
    showNA     = FALSE,      # drop NA legend entry,
    #textNA  = "Not Affected",
    legend.show = TRUE,
    fill.legend.labels = "Sentinel Sites",
    fill.legend.title = "",
    #colorNA = "white",
    #legend.is.portrait = FALSE
  ) +
  
  # bold ADM1 boundaries on top
  tm_shape(eth_adm1) +
  tm_borders(
    col       = "black",
    lwd       = 1
  ) +
  
  # layout adjustments
  tm_layout(
    title =  "Genomic Surveillance Sentinel sites in Ethiopia",
    title.size   = 1.2,
    title.position = c("center", "top"),
    #title.snap.to.legend = FALSE,
    #title.size = 1,
    legend.outside    = FALSE,
    legend.position = c(0.75, 0.75),
    #legend.outside.position  = c("right", "top"),
    # legend.title.size = 1.3,     # Smaller legend titles
    # legend.text.size = 1.1,      # Smaller legend text
    legend.xlab.size = 2,
    legend.ylab.size = 2,
    fill.legend.labels = "Sentinel Sites",
    fill.legend.title = "",
    #legend.show = FALSE,
    component.autoscale = FALSE,
    frame             = FALSE,
    bg.color          = "white",
    inner.margins = c(0.2, 0.15, 0.2, 0.4),  # enlarge right margin
    outer.margins = c(0, 0, 0, 0),
    fontfamily = "Arial"
  ) +
  
  # add scale bar & compass
  tm_scale_bar(
    position = c("left", "bottom"),
    text.size = 0.5
  ) +
  tm_compass(
    position = c("left", "top"),
    size     = 1
  )


tmap_save(senti_map, filename = "output/senti_map.png",
          width = 6, #Width in inches
          height = 6, #Height in inches
          dpi = 600 # DPI for high resolution
)


tmap_mode("plot")

senti_map2 <- 
  # ADM1 boundary as base
  tm_shape(eth_adm1) +
  tm_borders(col = "black", lwd = 1) +
  
  # ADM3 faint borders (optional, uncomment if needed)
  # tm_shape(eth_adm3) +
  # tm_borders(col = "grey90", lwd = 0.2) +
  
  # Sentinel sites polygons
  tm_shape(sent_site_sf) +
  tm_fill(
    col = "green",
    alpha = 0.7,
    legend.show = TRUE,
    legend.format = list(text = "Sentinel Sites")
  ) +
  
  # Optional: sentinel outline
  tm_borders(col = "darkgreen", lwd = 1) +
  tm_add_legend(
    type = "fill", 
    labels = "Sentinel Sites", 
    col = "green",
    border.col = "darkgreen"
  ) +
  
  # Layout settings
  tm_layout(
    title = "Genomic Surveillance Sentinel Sites in Ethiopia",
    title.size = 1.2,
    title.position = c("center", "top"),
    legend.show = TRUE,
    legend.title.size = 0,         # removes the legend title
    legend.text.size = 0.7,
    legend.position = c(0.75, 0.75),
    frame = FALSE,
    bg.color = "white",
    inner.margins = c(0.05, 0.1, 0.1, 0.1), 
    outer.margins = c(0,-0.005,0,-0.05),
    fontfamily = "Arial"
  ) +
  
  # Compass and scale bar
  tm_scale_bar(position = c("left", "bottom"), text.size = 0.5) +
  tm_compass(position = c("left", "top"), size = 1)

# Save the map
tmap_save(senti_map2, filename = "output/senti_map2.png",
          width = 6, height = 6, dpi = 600)





