label var num_families "Number of Families in 1993 (DANE)"

rename retornotierra land_return
label var land_return "Estimated land return follwing Salas (2014)"

rename aptitud land_aptitude
label var land_aptitude "Index of Land Suitability (IGAC)"

rename agua water
label var water "Index of Water availability (IGAC)"

label var drug_routes "Kms of Cocaine routes from Wright (2016)"

rename titulos_mineros_1990 mine_titles_90
label var mine_titles_90 "Number of Mining Titles, 1990 (Ministry of Energy)"

rename oil_prod_1998 oil_prod_98
label var oil_prod_98 "Oil production in 100 thousand barrels/day, 1998 (Ministry of Energy)"


rename funcionariosmunicipales1000hab local_bureaucracy_95
label var local_bureaucracy_95 "Number of local-level bureaucrats per 100, 1995 (Fundacion Social)"

rename gastomunicipalpercapita local_expenditure_95
label var local_expenditure_95 "Local public spending per 100, 1995 (Fundacion Social)"

rename notariaspc notaries_95
label var notaries_95 "Number of notaries per capita, 1995 (Fundacion Social)"

rename ofbancariaspc banks_95
label var banks_95 "Number of banks per capita, 1995 (Fundacion Social)"

rename ofcorpahorropc savings_95
label var savings_95 "Number of savings banks per capita, 1995 (Fundacion Social)"

rename centrossaludpc health_centers_95
label var health_centers_95 "Number of health centers per capita, 1995 (Fundacion Social)"

rename puestossaludpc health_posts_95
label var health_posts_95 "Number of notaries per capita, 1995 (Fundacion Social)"

rename planteleseducativospc schools_95
label var schools_95 "Number of schools per capita, 1995 (Fundacion Social)"

rename bibliotecaspblicaspc libraries_95
label var libraries_95 "Number of libraries per capita, 1995 (Fundacion Social)"

rename bomberospc fire_95
label var fire_95 "Number of fire stations per capita, 1995 (Fundacion Social)"

rename carcelespc jails_95
label var jails_95 "Number of jails per capita, 1995 (Fundacion Social)"

rename casasculturapc culture_95
label var culture_95 "Number of culture houses per capita, 1995 (Fundacion Social)"

rename ofinstrpublicospc instruments_95
label var instruments_95 "Number of public instruments offices per capita, 1995 (Fundacion Social)"

rename ofrecaudoimpuestospc tax_95
label var tax_95 "Number of tax collection offices per capita, 1995 (Fundacion Social)"



rename funcionariosestatales1000hab state_bureaucracy_95
label var state_bureaucracy_95 "Number of national-level bureaucrats per 100, 1995 (Fundacion Social)"

rename gastoestatalpercapita state_expenditure_95
label var state_expenditure_95 "National public spending per 100, 1995 (Fundacion Social)"

rename inspecpolpc police_95
label var police_95 "Police stations per capita, 1995 (Fundacion Social)"

rename puestospolpc police_posts_95
label var police_posts_95 "Police posts per capita, 1995 (Fundacion Social)"

rename juzgadospc courts_95
label var courts_95 "Courts per capita, 1995 (Fundacion Social)"

rename oftelecompc telecom_95
label var telecom_95 "Telecom offices per capita, 1995 (Fundacion Social)"

rename ofadpostalpc post_95
label var post_95 "Post offices capita, 1995 (Fundacion Social)"

rename ofcajaagrariapc agrarian_95
label var agrarian_95 "Agrarian bank offices per capita, 1995 (Fundacion Social)"

rename hospitalespc hospitals_95
label var hospitals_95 "Hospitals per capita, 1995 (Fundacion Social)"



rename stand_n_propietarios1_UAF owners_1_UAF
label var owners_1_UAF "Share of owners with 1 Family Agrarian Unit (CEDE)"

rename stand_n_propietarios2_UAF owners_2_UAF
label var owners_2_UAF "Share of owners with 2 Family Agrarian Units (CEDE)"

rename stand_n_propietarios3_UAF owners_3_UAF
label var owners_3_UAF "Share of owners with 3 Family Agrarian Units (CEDE)"


label var ratio_votes_dif_alc "Historical competitveness of mayoral elections (Pachon and Sanchez, 2014)"

label var ratio_votes_dif_asa "Historical competitveness of asembly elections (Pachon and Sanchez, 2014)"

label var ratio_votes_dif_pre "Historical competitveness of presidential elections (Pachon and Sanchez, 2014)"



rename acum_eventos_guerrilla guerrilla
label var guerrilla "Cumulative number of Guerrilla incidents (ODHVR)"


rename used_land_8 grazing
label var grazing "Hectares used for grazing (IGAC)"

rename cover_land_s9 grasses
label var grasses "Hectares of grasses (IGAC)"

rename cafe coffee
label var coffee "= 1 if coffee is grown (Ministry of Agriculture)"

rename ruggedness_index ruggedness
label var ruggedness "Land ruggedness index based on variation in elevation (IGAC)"

rename slope_mean slope
label var slope "Average terrain slope in degrees (IGAC)"

rename discapital_mean  dist_cap
label var dist_cap "Distance to department capital in Kms (IGAC)"

rename lnareaoficialhm2 area
label var area "Land area in square hectares (DANE)"

label var lat_mean "Latitude (IGAC)"

label var lon_mean "Longitude (IGAC)"

rename alt_mean elevation
label var elevation "Elevation in meters above sea level (IGAC)"

rename alt_mean2 elevation_2
label var elevation_2 "Square of elevation in meters above sea level (IGAC)"

rename rain_mean rainfall
label var rainfall "Rainfall in millimeters per year (IDEAM)"



label var landless_families "Proxy for the number of landless families (max(0, Families - Landowners)), (CEDE, IGAC, DANE)"

label var landowners_total "Number of Landowners in 2000, (CEDE, IGAC)"

label var landown_less_1 "Number of Landowners with up to 1 Ha in 2000, (CEDE, IGAC)"

label var landown_1_3 "Number of Landowners with 1-2 Ha in 2000, (CEDE, IGAC)"

label var landown_3_5 "Number of Landowners with 3-5 Ha in 2000, (CEDE, IGAC)"

label var landown_5_10 "Number of Landowners with 5-10 Ha in 2000, (CEDE, IGAC)"

label var landown_10_15 "Number of Landowners with 10-15 Ha in 2000, (CEDE, IGAC)"

label var landown_15_20 "Number of Landowners with 15-20 Ha in 2000, (CEDE, IGAC)"

label var landown_20_50 "Number of Landowners with 20-50 Ha in 2000, (CEDE, IGAC)"

label var landown_50_100 "Number of Landowners with 50-100 Ha in 2000, (CEDE, IGAC)"

label var landown_100_200 "Number of Landowners with 100-200 Ha in 2000, (CEDE, IGAC)"

label var landown_200_500 "Number of Landowners with 200-500 Ha in 2000, (CEDE, IGAC)"

label var landown_500_1000 "Number of Landowners with 500-1000 Ha in 2000, (CEDE, IGAC)"

label var landown_1000_2000 "Number of Landowners with 1000-2000 Ha in 2000, (CEDE, IGAC)"

label var landown_2000_plus "Number of Landowners with at least 2000 Ha in 2000, (CEDE, IGAC)"

label var convexity_L "Average Convexity of the Empirical Lorenz Curve (CEDE, IGAC)"

rename g_prop gini
label var gini "Land Gini among owners in 2000, (CEDE, IGAC)"




rename st_bov_nu2001 cattle_01
label var cattle_01 "Number of Cattle in 2001 (DANE)"

rename st_bov_nu2010 cattle_10
label var cattle_10 "Number of Cattle in 2010 (DANE)"

rename st_bov_nu2014 cattle_14
label var cattle_14 "Number of Cattle in 2014 (DANE)"

rename st_palm_s2002 palm_02
label var palm_02 "Hectares of palm planted in 2002 (Ministry of Agriculture)"

rename st_palm_s2011 palm_11
label var palm_11 "Hectares of palm planted in 2011 (Ministry of Agriculture)"




label var predicted_signal_strength_1 "Radio signal strength following Olken (2009), (Salcedo (2015))"


label var predicted_signal_strength_2 "Radio signal strength following Olken (2009) (alternative measure), (Salcedo (2015))"


label var free_space_signal "Topography-free radio signal strength following Olken (2009), (Salcedo (2015))"



label var educ_any_info "Fraction of adults with any education (DANE)"

rename educ_primary_orless educ_primary
label var educ_primary "Fraction of adults with Primary or less (DANE)"

rename educ_secundary_ormore educ_secundary
label var educ_secundary "Fraction of adults with Secondary or more (DANE)"















