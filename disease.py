import random

def transmit_hiv(G, transmission_probability_ftm,transmission_probability_mtf, current_step): #currently only hiv
    """
    For every partnership formed at `current_step`,
    attempt transmission from infected → susceptible.

    readability; person1 = person1, person2 = person2 whom they have a relationship with eachother
    """
    for person1, person2, attrs in G.edges(data=True):
        for a, b in ((person1, person2), (person2, person1)):
            if (G.nodes[a]['hiv_infection_status'] == 'I' and G.nodes[a]['gender'] == 'F' 
                and G.nodes[b]['hiv_infection_status'] == 'S'):
                if random.random() < transmission_probability_ftm:
                    G.nodes[b]['hiv_infection_status'] = 'I'
                    G.nodes[b]['hiv_infection_step'] = current_step
                    G.graph['hiv_total_infections'] += 1

            if (G.nodes[a]['hiv_infection_status'] == 'I' and G.nodes[a]['gender'] == 'M' 
                and G.nodes[b]['hiv_infection_status'] == 'S'):
                if random.random() < transmission_probability_mtf:
                    G.nodes[b]['hiv_infection_status'] = 'I'
                    G.nodes[b]['hiv_infection_step'] = current_step
                    G.graph['hiv_total_infections'] += 1




def transmit_flu(
    G,
    current_step,
    edge_beta=0.15,          # per-step flu transmission prob along a close/household-like edge
    hiv_multiplier=2.0,      # >1 increases risk for HIV+ susceptibles
    community_contacts=3,    # casual contacts per infectious person per step (well-mixed)
    community_beta=0.05      # per-contact prob for casual (non-edge) encounters
):
    """
    Influenza transmission for this time step.

    - Uses ALL active relationships (edges) each step, not just newly-formed ones.
    - Adds 'community mixing' so flu can jump outside the sexual-contact edges.
    - Susceptible people living with HIV ('hiv_infection_status' == 'I') have higher
      infection probability via `hiv_multiplier`.
    - Assumes node fields:
        G.nodes[i]['flu_status'] in {'S','E','I','R'}   (default 'S' if missing)
        G.nodes[i]['hiv_infection_status'] in {'S','I'} (existing in your model)
      On infection, we set:
        G.nodes[i]['flu_status'] = 'E'   (exposed; incubating)
        G.nodes[i]['flu_infection_step'] = current_step
        G.graph['flu_total_infections'] (counter; created if missing)
    """
    # --- helpers --------------------------------------------------------------
    def is_infectious(node_id):
        # Only 'I' shed flu; adjust if you later add presymptomatic infectiousness
        return G.nodes[node_id].get('flu_infection_status', 'S') == 'I'

    def is_susceptible(node_id):
        return G.nodes[node_id].get('flu_infection_status', 'S') == 'S'

    def hiv_multiplier_for(node_id):
        return hiv_multiplier if G.nodes[node_id].get('hiv_infection_status') == 'I' else 1.0

    def infect(node_id):
        G.nodes[node_id]['flu_infection_status'] = 'E'  # enter incubation; your progression will move E->I->R
        G.nodes[node_id]['flu_infection_step'] = current_step
        G.graph['flu_total_infections'] = G.graph.get('flu_total_infections', 0) + 1

    # To avoid multiple infections or double-processing, collect infections then apply once
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
    # Build quick lists for speed
    infectious = [n for n, d in G.nodes(data=True) if is_infectious(n)]
    susceptibles = [n for n, d in G.nodes(data=True) if is_susceptible(n)]

    # If no susceptibles remain, skip
    if susceptibles:
        for src in infectious:
            # Sample casual contacts without replacement, excluding self
            # (could bias by age, location, etc. later)
            pool = [x for x in susceptibles if x != src and x not in newly_infected]
            if not pool:
                continue
            k = min(community_contacts, len(pool))
            # random.sample requires k <= len(pool)
            contacts = random.sample(pool, k=k)
            for tgt in contacts:
                p = community_beta * hiv_multiplier_for(tgt)
                if random.random() < p:
                    newly_infected.add(tgt)

    # --- apply infections -----------------------------------------------------
    for n in newly_infected:
        infect(n)


def progress_flu(G, current_step, inc_period=1, inf_period=4):
    """
    E -> I after inc_period steps; I -> R after inf_period steps (from becoming I).
    """
    for n, attrs in G.nodes(data=True):
        status = attrs.get('flu_status', 'S')
        t_inf = attrs.get('flu_infection_step')

        if status == 'E' and t_inf is not None and current_step - t_inf >= inc_period:
            G.nodes[n]['flu_status'] = 'I'
            G.nodes[n]['flu_became_infectious_step'] = current_step

        elif status == 'I':
            t_I = attrs.get('flu_became_infectious_step', t_inf)
            if t_I is not None and current_step - t_I >= inf_period:
                G.nodes[n]['flu_status'] = 'R'
