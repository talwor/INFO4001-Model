import random

def transmit_hiv(G, transmission_probability_ftm, transmission_probability_mtf, current_step):
    new_today = 0
    for person1, person2, attrs in G.edges(data=True):
        for a, b in ((person1, person2), (person2, person1)):
            if G.nodes[a]['hiv_infection_status'] == 'I' and G.nodes[b]['hiv_infection_status'] == 'S':
                # base transmission prob
                p = (transmission_probability_ftm if G.nodes[a]['gender'] == 'F' else transmission_probability_mtf)

                #increase risk if either has influenza
                if (G.nodes[a]['flu_infection_status'] == 'I' or G.nodes[b]['flu_infection_status'] == 'I'):
                    p *= 2.0  #2 times the probability

                if random.random() < p:
                    G.nodes[b]['hiv_infection_status'] = 'I'
                    G.nodes[b]['hiv_infection_step'] = current_step
                    G.nodes[b]['hiv_ever_infected'] = True
                    new_today +=1
    return new_today




def transmit_flu(
    G,
    current_step,
    edge_beta,          # per-step flu transmission prob along a someone in a relationship
    hiv_multiplier,      # >1 increases risk for HIV+ susceptibles
    community_contacts,    # casual contacts per infectious person per step (well-mixed)
    community_beta  # per-contact prob for casual (non-edge) encounters. 5%
    ):     
    """
    Influenza transmission for this time step.

    - Uses ALL active relationships (edges) each step, not just newly-formed ones.
    - Adds 'community mixing' so flu can jump outside the sexual-contact edges.
    - Susceptible people living with HIV ('hiv_infection_status' == 'I') have higher
      infection probability via `hiv_multiplier`.

      On infection, we set:
        G.nodes[i]['flu_infection_status'] = 'E'   (exposed; incubating)
        G.nodes[i]['flu_infection_step'] = current_step
        G.graph['flu_total_infections'] (counter; created if missing)
    """
    # --- helpers --------------------------------------------------------------
    def is_infectious(node_id):
        # Only 'I' shed flu
        return G.nodes[node_id].get('flu_infection_status', 'S') == 'I'

    def is_susceptible(node_id):
        return G.nodes[node_id].get('flu_infection_status', 'S') == 'S'

    def hiv_multiplier_for(node_id):
        return hiv_multiplier if G.nodes[node_id].get('hiv_infection_status') == 'I' else 1.0

    def infect(node_id):
        G.nodes[node_id]['flu_infection_status'] = 'E'  # enter incubation; your progression will move E->I->R
        G.nodes[node_id]['flu_infection_step'] = current_step
        G.nodes[node_id]["flu_ever_infected"] = True

        G.graph['flu_total_infections'] = G.graph.get('flu_total_infections', 0) + 1

    # avoid multiple infections or double-processing, collect infections then apply once
    newly_infected = set()

    # --- 1) Transmission along edges (close/household-like contacts) ----------
    for person1, person2, attrs in G.edges(data=True):
        # two directions: person1->person2 and person2->person1
        if is_infectious(person1) and is_susceptible(person2) and person2 not in newly_infected:
            p = edge_beta * hiv_multiplier_for(person2)
            if random.random() < p:
                newly_infected.add(person2)

        if is_infectious(person2) and is_susceptible(person1) and person1 not in newly_infected:
            p = edge_beta * hiv_multiplier_for(person1)
            if random.random() < p:
                newly_infected.add(person1)

    # --- 2) Community mixing (casual contacts beyond the network edges) -------
    infectious = [n for n, d in G.nodes(data=True) if is_infectious(n)]
    susceptibles = [n for n, d in G.nodes(data=True) if is_susceptible(n)]

    # If no susceptibles remain, skip
    if susceptibles:  # only run if there are any susceptibles left
        for infector in infectious:
            # Build a pool of susceptible contacts (excluding the infector and those already infected this step)
            possible_contacts = [
                person for person in susceptibles
                if person != infector and person not in newly_infected
            ]
            if not possible_contacts:
                continue  # no one left to infect

            #number of casual encounters this infectious person will have this step
            num_contacts = min(community_contacts, len(possible_contacts))
            
            #randomly pick that many contacts from the pool
            chosen_contacts = random.sample(possible_contacts, k=num_contacts)

            for contact in chosen_contacts:
                #adjust transmission probability if contact is HIV-positive
                transmission_prob = community_beta * hiv_multiplier_for(contact)
                
                #infect with probability 'transmission_prob'
                if random.random() < transmission_prob:
                    newly_infected.add(contact)
                                


    # --- apply infections -----------------------------------------------------
    for person in newly_infected:
        infect(person)
    
    return len(newly_infected)  #cumulative incidence


def progress_flu(
    G, 
    current_step, 
    incubation_period=4, 
    infectious_period=7,
    immunity_days = 180
    ):
    """
    E is incubation, I is infected (E for exposed)
    E -> I after incubation_period steps; 
    I -> R after infectious_period steps (from becoming I).
    """
    for person, attrs in G.nodes(data=True):
        status = attrs.get('flu_infection_status', 'S')
        t_inf = attrs.get('flu_infection_step')

        if status == 'E' and t_inf is not None and current_step - t_inf >= incubation_period:
            G.nodes[person]['flu_infection_status'] = 'I'
            G.nodes[person]['flu_became_infectious_step'] = current_step

        elif status == 'I':
            t_I = attrs.get('flu_became_infectious_step', t_inf)
            if t_I is not None and (current_step - t_I) >= infectious_period:
                G.nodes[person]['flu_infection_status'] = 'R'
                G.nodes[person]['flu_recovered_step']   = current_step
                #heterogeneous waning (mean 180d, sd 30d; clamp to >=30)
                w = max(30, int(random.gauss(immunity_days, 30)))
                G.nodes[person]['flu_waning_days'] = w

        
        #changing from Recovered -> Susceptible 
        elif status == 'R':
            t_R = attrs.get('flu_recovered_step')
            if t_R is not None and current_step - t_R >= immunity_days:
                G.nodes[person]['flu_infection_status'] = 'S'
                #clear timestamps to allow reinfection
                G.nodes[person]['flu_infection_step'] = 0
                G.nodes[person]['flu_became_infectious_step'] = 0
                G.nodes[person]['flu_recovered_step'] = current_step

