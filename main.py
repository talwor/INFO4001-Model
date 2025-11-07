import numpy as np
import population_functions as pop
import relationship_functions as relationship
import networkx as nx
import disease_functions as disease
import random
import matplotlib.pyplot as plt


infected_hiv_counts = [20]
infected_flu_counts = [20]
new_flu_cases_per_day = [] 
new_hiv_cases_per_day = []

# cumulative_flu_counts = []


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

    # 3.simulate network growth over time
    total_steps = 730  # 730 days, 104 weeks, 2 years #at 365 for 1 year for nice graph
    for step in range(total_steps):
        print(step)
        disease.progress_flu(
            population,
            step,
            incubation_period=4,
            infectious_period=7
        )
        relationship.start_relationship(
            population,
            0.1429, #forimation_probablity #intercourse once a week
            0.7, #homophily
            step, #current step
            16, #min_age
            None
        )
        flu_new_today = disease.transmit_flu(
            population,
            step,
            edge_beta=0.15,
            hiv_multiplier=2.0,
            community_contacts=5,
            community_beta=0.1
        )
        new_flu_cases_per_day.append(flu_new_today)  

        hiv_new_today = disease.transmit_hiv(
            population,
            transmission_probability_mtf=1/1234, #based on research 1/1234
            transmission_probability_ftm=1/2380, #1/2380
            current_step=step
        ) 
        new_hiv_cases_per_day.append(hiv_new_today)  

        relationship.breakup(
            population,
            breakup_probability=0.02
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

    deg = dict(population.degree())
    print("max degree:", max(deg.values()))

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

    print("\ntotal num flu infections at end: "+ str(population.graph['flu_total_infections']))
    print("total flu recovered by the end: " + str(count_flu))

    print("===============")






    out = f"bourke_hiv_influenza_final_{step}d.graphml"

    nx.write_graphml(population, out)

    deg = dict(population.degree())
    nx.set_node_attributes(population, deg, "degree_now")



    weeks = list(range(0, total_steps+1))
    #HIV
    plt.figure(figsize=(8,4))
    plt.plot(weeks, infected_hiv_counts, linestyle='-')
    plt.xlabel("Day")
    plt.ylabel("Number of HIV-infected individuals")
    plt.title("HIV Prevalence over Time in Bourke Simulation")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    days = list(range(1, total_steps + 1))
    cum_hiv_cases = np.cumsum(new_hiv_cases_per_day)

    plt.figure(figsize=(9,4))
    plt.plot(days, cum_hiv_cases, label="HIV Cumulative cases (count)")
    plt.xlabel("Day")
    plt.ylabel("People ever infected")
    plt.title("Cumulative HIV incidence (Bourke model)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    #FLU
    plt.figure(figsize=(8,4))
    plt.plot(weeks, infected_flu_counts, linestyle='-')
    plt.xlabel("Day")
    plt.ylabel("Number of Flu-infected individuals (Prevalence)")
    plt.title("Influenza Prevalence over Time in Bourke Simulation")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    #build incidence & cumulative incidence vectors
    days = list(range(1, total_steps + 1))
    cum_flu_cases = np.cumsum(new_flu_cases_per_day)

    #plot daily new cases
    plt.figure(figsize=(9,4))
    plt.plot(days, new_flu_cases_per_day, linestyle='-')
    plt.xlabel("Day")
    plt.ylabel("New flu cases (incidence)")
    plt.title("Daily influenza incidence (Bourke model)")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    #plot cumulative incidence
    plt.figure(figsize=(9,4))
    plt.plot(days, cum_flu_cases, label="Flu Cumulative cases (count)")
    plt.xlabel("Day")
    plt.ylabel("People ever infected")
    plt.title("Cumulative influenza incidence (Bourke model)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
