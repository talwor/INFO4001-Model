import random
import networkx as nx

#functions to generate population

# def sample_from_distribution(distribution):
#     """
#     distribution: list of (value, probability) pairs summing to 1.
#     Returns one value sampled according to its probability.
#     """
#     r = random.random()
#     cumulative = 0.0
#     for value, prob in distribution:
#         cumulative += prob
#         if r < cumulative:
#             return value
#     return distribution[-1][0]

def generate_population(population_size, age_distribution, male_fraction, indigenous_fraction):
    """
    Builds a Graph where each node represents one person with attributes:
      - age
      - gender ('M' or 'F')
      - is_indigenous (True/False)
    """
    G = nx.Graph()
    G.graph['num_relationships_formed'] = 0
    G.graph['num_breakups'] = 0
    G.graph['total_infections'] = 0
    for person_id in range(population_size):
        age = sample_age(age_distribution)
        gender = 'M' if random.random() < male_fraction else 'F'
        is_indigenous = (random.random() < indigenous_fraction)
        G.add_node(person_id,
                   age=age,
                   gender=gender,
                   is_indigenous=is_indigenous,
                   infection_status="S", #S for susceptible, I for infected
                   infection_step=None) #infection_step = when the infection occurs
        
    adult_nodes = [
    node_id
    for node_id, attrs in G.nodes(data=True)
    if attrs['age'] >= 18
    ] #all people/nodes age 18 and over
        
    # e.g. choose 5 random seeds to start infected patient 0, patient 0 only can occur to an infectant over 18
    initial_infected_ids = random.sample(adult_nodes, 20)
    # mark them in the graph at t=0
    for seed in initial_infected_ids:
        G.nodes[seed]['infection_status'] = 'I'
        G.nodes[seed]['infection_step'] = 0

    return G


def sample_age(brackets):
    # pick a bracket first
    r = random.random()
    cum = 0
    for lo, hi, prob in brackets:
        cum += prob
        if r < cum:
            # then pick an actual age within [lo, hi]
            return random.randint(lo, hi)
    # fallback in case of rounding errors
    lo, hi, _ = brackets[-1]
    return random.randint(lo, hi)