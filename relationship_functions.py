import random

def find_eligible_partners(G, person_id, max_age_gap=5):
    """
    returns a list of node IDs that:
      - are not the same person
      - have opposite gender
      - are within max_age_gap years of age
    """
    person_age = G.nodes[person_id]['age']
    person_gender = G.nodes[person_id]['gender']
    eligible = []
    for candidate_id, attrs in G.nodes(data=True):
        if candidate_id == person_id:
            continue #not the same person allowed
        if attrs['gender'] == person_gender:
            continue #heterosexual mixing
        if abs(attrs['age'] - person_age) > max_age_gap:
            continue #enforce age gap limit
        eligible.append(candidate_id)
    return eligible

def start_relationship(G, formation_probability, homophily, current_step, transmission_probability,
                       min_age):
    """
    for each person, with probability formation_probability,
    attempt to form a new sexual partnership:
      - select from eligible partners
      - weight by whether they share indigenous status (homophily)
      - record edge attribute 'formed_step'
    """
    for person_id in list(G.nodes):
        if G.nodes[person_id]['age'] < min_age: #only 16 or older
            continue

        if G.degree(person_id) >= 1:  #monogamy
            continue

        if random.random() < formation_probability:
            # only consider candidates with no partners
            candidates = [
                p for p in find_eligible_partners(G, person_id)
                if G.degree(p) == 0
            ]
            if not candidates:
                continue

        #form relationship
        if random.random() < formation_probability:
            partners = find_eligible_partners(G, person_id)
            if not partners:
                continue
            #compute weights for homophily on indigenous status
            weights = []
            for partner_id in partners:
                same_status = (G.nodes[person_id]['is_indigenous'] ==
                               G.nodes[partner_id]['is_indigenous'])
                weights.append(homophily if same_status else (1 - homophily))
            total_weight = sum(weights)
            probabilities = [w / total_weight for w in weights]
            chosen_partner = random.choices(partners, probabilities)[0]

            if G.has_edge(person_id, chosen_partner):
                continue #prevent duplicate edges
            G.add_edge(person_id,chosen_partner,formed_step=current_step) #relationship formed
            G.graph['num_relationships_formed'] +=1

            #infecting, currently by chance.
            for a, b in ((person_id, chosen_partner), (chosen_partner, person_id)):
                if G.nodes[a]['infection_status'] == 'I' and G.nodes[b]['infection_status'] == 'S': 
                    if random.random() < transmission_probability:
                        G.nodes[b]['infection_status'] = 'I'
                        G.graph['total_infections'] +=1
                        G.nodes[b]['infection_step'] = current_step

           
def breakup(G, breakup_probability=0.01):
    """
    Each existing partnership has a chance to dissolve
    at each time step.
    """
    for u, v in list(G.edges()):
        if random.random() < breakup_probability:
            G.remove_edge(u, v)
            G.graph['num_breakups'] +=1
