#!/usr/bin/env python

#import .xml file 
import sys
from xml.dom import minidom
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np

path = sys.argv[1]
spim_root = ET.parse(path).getroot()

#Convert xml root to panda dataframe 
def xmltopandaDF(search_element='Spot',names = None):
    if names is None:
        for spot in spim_root.iter(search_element):
            names = spot.keys()
            break
    
    df = pd.DataFrame(columns = names)
    #Find particular elements of interest and iterate trough list of all subelements under the root
    for i,spots in enumerate(spim_root.iter(search_element)):
        #fill panda df row by row 
        df.loc[i] = spots.attrib
    return df

spots_df = xmltopandaDF("Spot")
edges_df = xmltopandaDF("Edge")
spot_df_outfile = ".".join(path.split(".")[:(len(path.split("."))-1)] + ["spot_df.csv"])
edge_df_outfile = ".".join(path.split(".")[:(len(path.split("."))-1)] + ["edge_df.csv"])
spots_df.to_csv(spot_df_outfile)
edges_df.to_csv(edge_df_outfile)

