# -*- coding: utf-8 -*-
"""
Created in Dec 2020
@author: Tom Nelson
"""

import json
MONOMER_DB = "monomer_org.json"
DEFAULT_POLYMER_TYPE = "PEPTIDE"

class Monomer:
    '''Class for representing monomers'''
    def __init__(self, symbol, polymerType=DEFAULT_POLYMER_TYPE):
        with open(MONOMER_DB, "r") as dbfile:
            data = json.load(dbfile)
        found = False
        for item in data:
            if item['symbol'] == symbol and item['polymerType'] == polymerType:
                found = True
                self.author = item['author']
                self.createDate = item['createDate']
                self.id = item['id']
                self.molfile = item['molfile']
                self.monomerType = item['monomerType']
                self.name = item['name']
                self.naturalAnalog = item['naturalAnalog']
                self.polymerType = item['polymerType']
                self.rgroups = item['rgroups']
                self.smiles = item['smiles']
                self.symbol = item['symbol']
        if not found:
            raise NameError

class SimplePolymer:
    
    def __init__(self, monomerList, polymerType):
        self.monomerList = monomerList
        self.polymerType = polymerType
        self.simplePolymer = []
        for item in self.monomerList:
            self.simplePolymer.append(Monomer(item, self.polymerType))



### Testing ###
myMonomer1 = Monomer('C', "RNA")
myMonomer1.smiles
myMonomer1.rgroups
myMonomer1.name

myMonomer2 = Monomer('C', "PEPTIDE")
myMonomer2.smiles
myMonomer2.rgroups
myMonomer2.name

myMonomer3 = Monomer('Aib')
myMonomer3.name
myMonomer3.naturalAnalog

myMonomer4 = Monomer('foobar')

spam = SimplePolymer(['M','A','L'], "PEPTIDE")
spam.simplePolymer[0].symbol
spam.simplePolymer[0].name
spam.simplePolymer[0].smiles


