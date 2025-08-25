import population_functions as pop
import relationship_functions as relationship
import networkx as nx
import disease
import random

import matplotlib.pyplot as plt

infected_hiv_counts = [20]
infected_flu_counts = [20]


"""
Simulation in code below
"""



if __name__ == "__main__":
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
    nx.write_graphml(population, "bourke_hiv_influenza_start.graphml")


    #adding how many people are infected w/ HIV
    count = 0 
    for person_id, attrs in population.nodes(data=True):
        if attrs["hiv_infection_status"] == "I":
            count +=1
    initialisation_count_hiv = count
    
    count = 0 
    for person_id, attrs in population.nodes(data=True):
        if attrs["flu_infection_status"] == "I":
            count +=1
    initialisation_count_flu = count


    print("numbreakups at initialization (should be 0): "+ str(population.graph['num_breakups']))
    print("num relationships at initialization (should be 0): "+ str(population.graph['num_relationships_formed']))




    # 3. Simulate network growth over time
    total_steps = 730  # 730 days, 104 weeks, 2 years
    for step in range(total_steps):
        relationship.start_relationship(
            population,
            0.01, #formation_probability
            0.7, #homophily
            step, #current step
            16 #min_age
        )
        disease.progress_flu(
        population,
        step,
        incubation_period=4,
        infectious_period=7
        )
        disease.transmit_flu(
            population,
            step,
            edge_beta=0.15,
            hiv_multiplier=2.0,
            community_contacts=5,
            community_beta=0.05
        )
        disease.transmit_hiv(
            population,
            transmission_probability_mtf=1/1234, #based on research 1/1234
            transmission_probability_ftm=1/2380, #1/2380
            current_step=step
        ) 
        relationship.breakup(
            population,
            breakup_probability=0.01 #0.001
        )
        pop.apply_recovery(
            population,
            step
        )
        
        #matplot lib
        count_hiv = sum(1 for _, attrs in population.nodes(data=True)
        if attrs["hiv_infection_status"] == "I")
        infected_hiv_counts.append(count_hiv)

        count_flu = sum(1 for _, attrs in population.nodes(data=True)
        if attrs["flu_infection_status"] == "I")
        infected_flu_counts.append(count_flu)

    count_hiv = 0 
    for person_id, attrs in population.nodes(data=True):
        if attrs["hiv_infection_status"] == "I":
            count_hiv +=1

    count_flu = 0 
    for person_id, attrs in population.nodes(data=True):
        if attrs["flu_infection_status"] == "I":
            count_flu +=1



    """
    initialisation stats:
    """
    print("\n====INITIALISATION====")
    print("hiv infected count at initialization: " + str(initialisation_count_hiv))
    print("flu infected count at initialization: " + str(initialisation_count_flu))


    """
    runtime
    """
    print("\n====RUNTIME====")
    print("total amt of infected with HIV: " + str(count_hiv))
    print("numbreakups: "+ str(population.graph['num_breakups']))
    print("num hiv relationships formed in total over "+ str(total_steps)+ " weeks: "+ str(population.graph['num_relationships_formed']))

    print("total amt of infected with flu: " + str(count_flu))

    """
    end
    """

    count_hiv = 0 
    for person_id, attrs in population.nodes(data=True):
        if attrs["hiv_infection_status"] == "R":
            count_hiv +=1

    count_flu = 0 
    for person_id, attrs in population.nodes(data=True):
        if attrs["flu_infection_status"] == "R":
            count_flu +=1


    print("\n====END====")
    print("total num hiv infections at end: "+ str(population.graph['hiv_total_infections']))
    print("total hiv recovered by the end: " + str(count_hiv))
    print("nodes:", population.number_of_nodes(), "edges:", population.number_of_edges())
    print(list(population.edges(data=True))[:5]) 

    print("\ntotal num flu infections at end: "+ str(population.graph['flu_total_infections']))
    print("total flu recovered by the end: " + str(count_flu))

    nx.write_graphml(population, "bourke_hiv_influenza_final.graphml")



    weeks = list(range(0, total_steps+1))
    #HIV
    plt.figure(figsize=(8,4))
    plt.plot(weeks, infected_hiv_counts, marker='o')
    plt.xlabel("Day")
    plt.ylabel("Number of HIV-infected individuals")
    plt.title("HIV Prevalence over Time in Bourke Simulation")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    #FLU
    plt.figure(figsize=(8,4))
    plt.plot(weeks, infected_flu_counts, marker='o')
    plt.xlabel("Day")
    plt.ylabel("Number of Flu-infected individuals (I)")
    plt.title("Influenza Prevalence over Time in Bourke Simulation")
    plt.grid(True)
    plt.tight_layout()
    plt.show()



    # now `population` holds  dynamic sexual-contact network,
    # with edges stamped by 'formed_step' when they appeared.
