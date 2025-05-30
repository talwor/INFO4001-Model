import population_functions as pop
import relationship_functions as relationship
import networkx as nx
import random


# 1. Fill these from ABS data:
population_size = 2340
age_distribution = [
    (0, 4,   178/population_size),
    (5, 9,   172/population_size),
    (10, 14, 165/population_size),
    (15, 19, 134/population_size),
    (20, 24, 134/population_size),
    (25, 29, 150/population_size),
    (30, 34, 168/population_size),
    (35, 39, 148/population_size),
    (40, 44, 122/population_size),
    (45, 49, 128/population_size),
    (50, 54, 173/population_size),
    (55, 59, 164/population_size),
    (60, 64, 158/population_size),
    (65, 69, 102/population_size),
    (70, 74, 110/population_size),
    (75, 79,  67/population_size),
    (80, 84,  46/population_size),
    (85, 120, 34/population_size),  #assume 85–120 covers 85+
]
male_fraction        = 1164 / population_size
indigenous_fraction  = 708 / population_size

# 2. Generate the population graph
population = pop.generate_population(
    population_size,
    age_distribution,
    male_fraction,
    indigenous_fraction
)



count = 0 
for person_id, attrs in population.nodes(data=True):
    if attrs["infection_status"] == "I":
        count +=1
print("infected count at initialization: " + str(count))
print("numbreakups at initialization (should be 0): "+ str(population.graph['num_breakups']))
print("num relationships at initialization (should be 0): "+ str(population.graph['num_relationships_formed']))





# 3. Simulate network growth over time
total_steps = 120  # e.g. months
for step in range(total_steps):
    relationship.start_relationship(
        population,
        0.05, #formation_probability
        0.7, #homophily
        step, #current step
        0.1, #transmission_probability
        16 #min_age
    ) 
    relationship.breakup(
        population,
        breakup_probability=0.01
    )

count = 0 
for person_id, attrs in population.nodes(data=True):
    if attrs["infection_status"] == "I":
        count +=1
        
print("total amt of infected: " + str(count))
print("numbreakups: "+ str(population.graph['num_breakups']))
print("num relationships formed in total over 120 steps: "+ str(population.graph['num_relationships_formed']))

print("total num infections at end: "+ str(population.graph['total_infections']))


# now `population` holds  dynamic sexual-contact network,
# with edges stamped by 'formed_step' when they appeared.
