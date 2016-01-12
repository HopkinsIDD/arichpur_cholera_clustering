######### arichpur_cholera_clustering

######### q1q4iccctrl = household level exposures from spatially matched households
# "dataid" = household id
# "CODE" = matched-sest id
# "X","Y","longitude","latitude" = household GPS coordinates blocked for privacy concerns
# household level exposures ("use of supplied water","use of tubewell water","distance to the\n closest water source", 
# "intermittent water supply","boiling water","handsoap availability",
# "type of a latrine","household density","sharing a latrine") 
# correspond to the variables ("q11_supply","q11_tube","dist_more10", 
# "q4water2","boilnewst", "soap",
# "defecate_pitl","pplrm_di","defsharing_only")


######## q2icc_nohh = individual level exposures from spatially matched households
# "datamemberid" = participant id
# "dataid" = household id
# "CODE" = matched-sest id
# "X","Y","longitude","latitude" = household GPS coordinates blocked for privacy concerns
# individual level exposures ("Feeding a child with hand in the past week",
# "Eating meals prepared over 2 hours before consumption in the past week",
# "Drinking water outside the home in the past week",
# "Eating fresh cut fruit or vegetables outside the home in the past week",
# "Drinking tea outside the home in the past week") 
# correspond to the variables ("handfeed_never","prep2hrs_never","q40_bc_never","q41_4_bc_never","q41_5_bc_never")

######## arichpur_clustering_sourcecodes
# source codes for generating within matched-sets concordance and between matched-sets concordance 
# of both household and individual level exposures

######## co-occurrence_sourcecodes
# source codes for generating co-occurrence of exposures 
# of both household and individual level exposures

######## sample codes
# examples for running source code files

######## ICC.do
# codes for generating ICC values for both household and individual level exposures
