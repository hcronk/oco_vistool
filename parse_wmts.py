#
# the input file, 'wmts.xml', is the "GIBS Capabilities Document", 
# can be accessed with a browser at this URL:
#
# https://gibs.earthdata.nasa.gov/wmts/epsg4326/best/wmts.cgi?SERVICE=WMTS&REQUEST=GetCapabilities
#
# The returned XML document should be saved to disk and processed with this script.
#

import xml.etree.ElementTree as ET
import pandas as pd

tree = ET.parse("wmts.xml")
root = tree.getroot()
layer_list = list()
filename_list = list()

for contents in root.findall('{http://www.opengis.net/wmts/1.0}Contents'):
    for layer in contents.findall('{http://www.opengis.net/wmts/1.0}Layer'):
        for identifier in layer.findall('{http://www.opengis.net/ows/1.1}Identifier'):
            layer_list.append(identifier.text.strip())
        filename = layer.find('{http://www.opengis.net/ows/1.1}Metadata')
        if (filename is not None):
            filename_list.append(filename.attrib.get('{http://www.w3.org/1999/xlink}href'))
        else:
            filename_list.append('Null')

layer_xml_df = pd.DataFrame(list(zip(layer_list, filename_list)), columns = ['Name', 'Url'])
layer_xml_df.to_csv('Layer_xml.csv', index = False)
