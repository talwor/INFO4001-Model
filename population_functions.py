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
      etc.
    """
    G = nx.Graph()
    G.graph["num_relationships_formed"] = 0
    G.graph["num_breakups"] = 0
    G.graph["hiv_total_infections"] = 0
    G.graph["flu_total_infections"] = 0
    for person_id in range(population_size):
        age = sample_age(age_distribution)
        gender = 'M' if random.random() < male_fraction else 'F'
        is_indigenous = (random.random() < indigenous_fraction)
        G.add_node(person_id,
                   age=age,
                   gender=gender,
                   is_indigenous=is_indigenous,
                   hiv_infection_status="S", #S for susceptible, I for infected, R for recovered
                   hiv_infection_step=0, #infection_step = when the infection occurs
                   hiv_ever_infected=False,
                   flu_infection_status="S",
                   flu_infection_step=0,
                   flu_recovered_step=0,
                   flu_ever_infected=False) 
        
    adult_nodes = [
    node_id
    for node_id, attrs in G.nodes(data=True)
    if attrs['age'] >= 18
    ] #all people/nodes age 18 and over

    # e.g. choose 20 random seeds to start infected patient 0, patient 0 only can occur to an infectant over 18
    initial_infected_ids = random.sample(adult_nodes, 20)
    # mark them in the graph at t=0
    for seed in initial_infected_ids:
        G.nodes[seed]['hiv_infection_status'] = 'I'
        G.nodes[seed]['hiv_infection_step'] = 0 

    flu_seed_count = 20   # choose number of flu patient-zeros
    flu_seed_ids = random.sample(list(G.nodes), flu_seed_count)
    for seed in flu_seed_ids:
        G.nodes[seed]['flu_infection_status'] = 'I'             # start infectious
        G.nodes[seed]['flu_infection_step'] = 0
        G.nodes[seed]['flu_ever_infected'] = True        # <-- add
        G.graph['flu_total_infections'] += 1             # <-- add
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


# constants
RECOVERY_DAYS = 180  #change to ceil(60/7) if timestep is a WEEK

def apply_recovery(G, current_day, status_key="hiv_infection_status",
                   inf_step_key="hiv_infection_step", rec_step_key="hiv_recovery_step",
                   recovery_days=RECOVERY_DAYS):
    """
    Turn I -> R if (current_day - infection_step) >= recovery_days.
    """
    for n, attrs in G.nodes(data=True):
        if attrs[status_key] == "I" and (current_day - attrs[inf_step_key]) >= recovery_days:
            attrs[status_key] = "R"
            attrs[rec_step_key] = current_day
