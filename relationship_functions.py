import random

def find_eligible_partners(G, person_id,max_degree=10, max_age_gap=10):
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
        if candidate_id == person_id: #not self
            continue
        if attrs['gender'] == person_gender: #heterosexual
            continue
        if abs(attrs['age'] - person_age) > max_age_gap or attrs['age'] < 16: 
            continue
        # allow multiple partners; cap if max_degree is set
        if max_degree is not None and G.degree(candidate_id) >= max_degree:
            continue
        eligible.append(candidate_id)
    return eligible

def start_relationship(G, formation_probability, homophily, current_step, min_age, max_degree):
    """
    for each person, with probability formation_probability,
    attempt to form a new sexual partnership:
      - select from eligible partners
      - weight by whether they share indigenous status (homophily)
      - record edge attribute 'formed_step'
    """
    for person_id in list(G.nodes):
        if G.nodes[person_id]['age'] < min_age:
            continue
        #allow multiple partners; cap if max_degree is set
        if max_degree is not None and G.degree(person_id) >= max_degree:
            continue

        if random.random() < formation_probability:
            partners = find_eligible_partners(G, person_id, max_age_gap=10, max_degree=max_degree)
            if not partners:
                continue
            weights = []
            for partner_id in partners:
                same = (G.nodes[person_id]['is_indigenous'] == G.nodes[partner_id]['is_indigenous'])
                weights.append(homophily if same else (1 - homophily))
            total = sum(weights)
            probs = [w/total for w in weights]
            chosen = random.choices(partners, probs)[0]
            G.add_edge(person_id, chosen, formed_step=current_step)
            G.graph['num_relationships_formed'] += 1
           
def breakup(G, breakup_probability):
    """
    each existing partnership has a chance to dissolve
    at each time step. 
    """
    for person1, person2 in list(G.edges()):
        if random.random() < breakup_probability:
            G.remove_edge(person1, person2)
            G.graph['num_breakups'] +=1
