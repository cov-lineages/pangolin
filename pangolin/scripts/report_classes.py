#!/usr/bin/env python3
import datetime as dt
from collections import Counter

class taxon():
    
    def __init__(self, record_id, lineage):
    
        self.id = record_id
        self.lineage = lineage
        
        self.country = self.id.split("|")[1]

        self.get_date_loc()
        
        
    def get_date_loc(self):
        
        self.date = self.id.split("|")[2]
        
        date_bits = self.date.split("-")

        if len(date_bits) == 3:
            self.date_dt = dt.date(int(date_bits[0]), int(date_bits[1]), int(date_bits[2]))
        else:
            self.date_dt = "NA"

class lineage():
    
    def __init__(self, name, taxa):
        
        self.id = name
        
        if self.id == "":
            self.new = True
        else:
            self.new = False

        self.taxa = taxa
        
        self.dates = []
        #self.epiweeks = []
        
        self.locations = []
        self.main_locs = []
        
        self.country_freqs = {}
        
        self.get_date_loc_info()
        self.get_most_common_country()
        
        

        
    def get_date_loc_info(self):
        
        current_day = dt.date.today()
        
        for tax in self.taxa:
            if tax.date_dt != "NA":
                self.dates.append(tax.date_dt)
                self.locations.append(tax.country)

        self.date_counts = Counter(self.dates)
        self.loc_counts = Counter(self.locations)
                
        if self.dates == []:
            #print(self.id)
            pass
        else:   
            self.mrd = max(self.dates)
            self.oldest = min(self.dates)
            
            self.pretty_mrd = self.mrd.strftime('%B-%d')
            self.pretty_oldest = self.oldest.strftime('%B-%d')

            self.length = (self.mrd - self.oldest).days
            
        self.last_sampled = (current_day - self.mrd).days

        
    def get_most_common_country(self):
        
        total = len(self.taxa)
        
        for country, count in self.loc_counts.items():
            
            frequency = count/total
            
            self.country_freqs[country] = count/total
            
            
        all_3 = self.loc_counts.most_common(3)
        self.main_locs = [i[0] for i in all_3]
