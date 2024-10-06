import csv
import pandas as pd
import numpy as np
import json

# Define the ExoPlanet class
#  	Name	Type	Detection Method	Mass
class ExoPlanet:
    def __init__(self, name, type, detection_method, mass, radius, flux, Tsurf, period, edistance, age, ESI, found):
        self.name = name
        self.type = type
        self.detection_method = detection_method
        self.mass = mass
        self.radius = radius
        self.flux = flux
        self.Tsurf = Tsurf
        self.period = period
        self.edistance = edistance
        self.age = age
        self.ESI = ESI
        self.found = found


# Define a hash table where the key is the solar system name and value is a list of ExoPlanet objects
Habitable_List = {}

# Read the CSV file and populate the hash table
with open('hwc_table_all.csv', 'r') as file:
    reader = csv.reader(file)
    header = next(reader)  # Skip the header row if present
    
    for row in reader:
        # Assuming columns: name, star_name, mass, radius, distance, discovery_year
        name = row[0]

        #need to input correct colum vals cos these are not it
        type = row[1]
        detection_method = row[2]
        mass = float(row[3])
        radius = float(row[4])
        flux = float(row[5])
        Tsurf = row[6]
        period = float(row[7])
        edistance = float(row[8])
        age = row[9]
        ESI = float(row[10])
        found = 0
        
        
        # Create an ExoPlanet object
        current = ExoPlanet(name, type, detection_method, mass, radius, flux, Tsurf, period, edistance, age, ESI, found)
        
        # Store the exoplanet in the Exo_List hash table using the name as the key
        if name not in Habitable_List:
            Habitable_List[name] = current  # Initialize a new list for planets in this solar system

# Define the constants in SI units BE CONSISTENT USING SI UNITS
# Gravitational constant
G = 6.67*10**(-11)
# Solar mass
M = 1.989*10**(30)
# Solar luminosity
L = 3.826*10**(26)
# Solar radius, in earth masses
r = 6.9599*(10**8) / 6378137
# Boltzmann constant (the little k used in the distance fo)
k = 1.380649*10**(-23)
# Define the ExoPlanet class

# Spectral Class Masses
mass_G = M
mass_K = M * 0.746
mass_M = M * 0.442

# Spectral type radii
rad_G = r
rad_K = r * 0.9
rad_M = r * 0.5

#  	Name	Type	Detection Method	Mass
class SolarSystem:
    def __init__(self, starName, radius, starMass, allExo, numPlanets):
        self.name = starName
        self.radius = radius
        self.starMass = starMass
        self.allExo = allExo
        self.numPlanets=numPlanets

    def __repr__(self):
        return (f"SolarSystem(name={self.name}, radius={self.radius}, starMass={self.starMass}, "
                f"allExo={self.allExo}, numPlanets={self.numPlanets})")

class Solar_Sys_Expo:
    def __init__(self, name, radius, period, home_star, star_temp_eff, dist_from_star, orb_speed):
        self.name=name
        self.radius = radius
        self.period = period
        self.home_star = home_star
        self.star_temp_eff=star_temp_eff
        self.dist_from_star=dist_from_star
        self.orb_speed = orb_speed
        self.is_habitable = 0
        
    def __repr__(self):
        return (f"Solar_Sys_Expo(name={self.name}, radius={self.radius}, "
                f"period={self.period}, home_star={self.home_star},dist_from_star={self.dist_from_star})")
    

# Read the CSV file and replace any all-empty columns with 0
df = pd.read_csv('calculations.csv').fillna(0)

Exo_List = []

# Perform the groupby operation and calculate the mean
merged_df = df.groupby(['pl_name', 'hostname']).agg(
    pl_rade=('pl_rade', lambda x: np.nanmean(x)), 
    pl_masse=('pl_masse', lambda x: np.nanmean(x)),
    pl_orbper=('pl_orbper', lambda x: np.nanmean(x)),
    st_teff=('st_teff', lambda x: np.nanmean(x))
).reset_index()


for index, row in merged_df.iterrows():
    data = row.to_dict()  # Convert the row to a dictionary
    pl_name = data['pl_name']
    hostname = data['hostname']
    pl_rade = data['pl_rade']
    pl_masse = data['pl_masse']
    pl_orbper = data['pl_orbper']
    st_teff = data['st_teff']
    Exo_List.append(Solar_Sys_Expo(pl_name, pl_rade, pl_orbper, hostname, st_teff, 0, 0))

allSolarSystems = {}
# Iterate through list of exoplanets, putting exoplanets that have the necessary data 
# into a hash table for a particular solar system. Delete if null period, delete when moved
# AND calculating the distance from the star for each planet based on period
for Item in Exo_List:
    if Item.period == 0:
        if Item.name in Habitable_List:
            del Habitable_List[Item.name]
        continue
    else:
        # If the solar system hasn't yet been pushed into hashtable
        if Item.name in Habitable_List:
            Item.is_habitable = 1
            Habitable_List[name].found=1
        if Item.home_star not in allSolarSystems:
            starName= Item.home_star
            # Set radius and mass based on temperature (spectral type)
            if Item.star_temp_eff <= 3750:
                radius = rad_M
                starMass = mass_M
            elif Item.star_temp_eff <= 5240:
                radius = rad_K
                starMass = mass_K
            else:
                radius = rad_G
                starMass = mass_G
            Item.dist_from_star = (G * starMass * (Item.period**2) / (4 * np.pi**2)) ** (1/3)
            Item.orb_speed = 2 * np.pi * Item.dist_from_star / (Item.period*86400)
            allExo = [Item] 
            numPlanets=1
            tempSolar=SolarSystem(starName, radius, starMass, allExo, numPlanets)
            allSolarSystems[starName]= tempSolar
            
        else:
            currentStarMass=allSolarSystems[starName].starMass
            #use mass to calculate dist_from_star here in m
            Item.dist_from_star = (G * currentStarMass * Item.period**2 / (4 * np.pi**2)) ** (1/3)
            #use (something) to calculate orbital speed
            Item.orb_speed = 2 * np.pi * Item.dist_from_star / (Item.period*86400)
            allSolarSystems[starName].allExo.append(Item)
            allSolarSystems[starName].numPlanets+=1


keys = list(Habitable_List.keys())  # Create a copy of the keys
i = 0  # Initialize index for while loop




def solar_system_to_dict(solar_system):
    return {
        'name': solar_system.name,
        'radius': solar_system.radius,
        'starMass': solar_system.starMass,
        'allExo': [solar_sys_expo_to_dict(exo) for exo in solar_system.allExo],
        'numPlanets': solar_system.numPlanets
    }

def solar_sys_expo_to_dict(exo):
    return {
        'name': exo.name,
        'radius': exo.radius,
        'period': exo.period,
        'home_star': exo.home_star,
        'dist_from_star': exo.dist_from_star
    }


def habitable_exo_to_dict(exo):
    return {
        'name': exo.name,
        'type': exo.type,
        'detection_method': exo.detection_method,
        'mass': exo.mass,
        'radius': exo.radius,
        'flux': exo.flux,
        'Tsurf': exo.Tsurf,
        'period': exo.period,
        'edistance': exo.edistance,
        'age': exo.age,
        'ESI': exo.ESI,
        'found': exo.found
    }

solar_systems_dict = {key: solar_system_to_dict(value) for key, value in allSolarSystems.items()}




exoplanet_names = []

for star in solar_systems_dict.values():
    for exoplanet in star["allExo"]:
        exoplanet_names.append(exoplanet["name"])

while i < len(keys):  # Iterate over the list of keys
    key = keys[i]
    if key not in exoplanet_names:
        del Habitable_List[key]  # Delete key from the original dictionary
    i += 1  # Move to the next key




habitable_worlds_dict = {key: habitable_exo_to_dict(value) for key, value in Habitable_List.items()}

# Convert to JSON
json_output = json.dumps(solar_systems_dict, indent=4)
with open('solar_systems.json', 'w') as json_file:
    json_file.write(json_output)


json_output = json.dumps(habitable_worlds_dict, indent=4)
with open('habitable_worlds.json', 'w') as json_file:
    json_file.write(json_output)
