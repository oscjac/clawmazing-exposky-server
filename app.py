import csv

# Define the constants in SI units BE CONSISTENT USING SI UNITS
# Gravitational constant
G = 6.67*10**(-11)
# Solar mass
M = 1.989*10**(30)
# Solar luminosity
L = 3.826*10**(26)
# Boltzmann constant (the little k used in the distance fo)
k = 1.380649*10**(-23)
# Define the ExoPlanet class
#  	Name	Type	Detection Method	Mass
class SolarSystem:
    def __init__(self, starName, radius, allExo, numPlanets):
        self.name = starName
        self.radius = radius
        self.allExo = allExo
        self.numPlanets=numPlanets

    def __repr__(self):
        return (f"SolarSystem(name={self.name}, radius={self.radius}, "
                f"allExo={self.allExo}, numPlanets={self.numPlanets})")


class Solar_Sys_Expo:
    def __init__(self, name, radius, period, home_star, dist_from_star):
        self.name=name
        self.radius = radius
        self.period = period
        self.dist_from_star = dist_from_star
        
    def __repr__(self):
        return (f"Solar_Sys_Expo(name={self.name}, radius={self.radius}, "
                f"period={self.period}, home_star={},dist_from_star={self.dist_from_star})")



# Read the CSV file and populate the hash table
with open('change to file name', 'r') as file:
    reader = csv.reader(file)
    header = next(reader)  # Skip the header row if present

    # Define a hash table where the key is the exoplanet name and it maps to its data
    Exo_List = {}
    
    
    for row in reader:
        # Assuming columns: name, star_name, mass, radius, distance, discovery_year
        name = row[0]

        #need to input correct colum vals cos these are not it
        radius = float(row[1])
        period = float(row[3])
        dist_from_star = 20
        
        
        # Create an ExoPlanet object
        current = ExoPlanet(name, radius, period, dist_from_star)
        
        # Store the exoplanet in the Exo_List hash table using the name as the key
        if name not in Exo_List:
            Exo_List[name] = current  # Initialize a new list for planets in this solar system
        else :
            for i in range(2):
            if nnot Exo_List[name][i+1]:

# Iterate through list of exoplanets, putting exoplanets that have the necessary data 
# into a hash table for a particular solar system. Delete if null period, delete when moved
# AND calculating the distance from the star for each planet based on period
while Exo_List:
    

    

        
        

