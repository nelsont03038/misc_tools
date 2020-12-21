# -*- coding: utf-8 -*-
"""
Created in Dec 2020

@author: Tom Nelson
"""

import json
with open("monomer_org.json", "r") as read_file:
    data = json.load(read_file)

def get_info_from_symbol(symbol, polymerType='PEPTIDE'):
	for item in data:
		if item['symbol'] == symbol and item['polymerType'] == polymerType:
			return item['symbol'], item['naturalAnalog'], item['name']

get_info_from_symbol('C')
get_info_from_symbol('C', "RNA")
get_info_from_symbol('Aib')
get_info_from_symbol('Dab')




