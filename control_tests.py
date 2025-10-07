
import argparse
import random
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

import population_functions as pop
import relationship_functions as relationship
import disease



"""
HIV only, 2 years (730 days)
python control_tests.py --mode hiv_only --days 730 --seed 42

Flu only
python control_tests.py --mode flu_only --days 730

Both (baseline you’ve been using)
python control_tests.py --mode both --days 730

"""
# ----------------------------
# Utilities to (re)initialise
# ----------------------------
def reset_flu(G):
    """Clear all flu state from the graph and counters so we can run an HIV-only scenario."""
    for n, d in G.nodes(data=True):
        d['flu_infection_status'] = 'S'
        d['flu_infection_step'] = 0
        d['flu_recovered_step'] = 0
        d['flu_ever_infected'] = False
        d.pop('flu_became_infectious_step', None)
        d.pop('flu_waning_days', None)
    G.graph['flu_total_infections'] = 0


def reset_hiv(G):
    """Clear all HIV state from the graph and counters so we can run a Flu-only scenario."""
    for n, d in G.nodes(data=True):
        d['hiv_infection_status'] = 'S'
        d['hiv_infection_step'] = 0
        d['hiv_ever_infected'] = False
        d.pop('hiv_recovery_step', None)
    G.graph['hiv_total_infections'] = 0


def seed_only_hiv(G, num_seeds=20, min_age=18):
    """Seed HIV in adults only."""
    adult_nodes = [n for n, a in G.nodes(data=True) if a.get('age', 0) >= min_age]
    seeds = random.sample(adult_nodes, min(num_seeds, len(adult_nodes)))
    for s in seeds:
        G.nodes[s]['hiv_infection_status'] = 'I'
        G.nodes[s]['hiv_infection_step'] = 0
        G.nodes[s]['hiv_ever_infected'] = True
    # counter not incremented on seeding by design; track infections during runtime


def seed_only_flu(G, num_seeds=20):
    """Seed flu (can be any ages)."""
    nodes = list(G.nodes)
    seeds = random.sample(nodes, min(num_seeds, len(nodes)))
    for s in seeds:
        G.nodes[s]['flu_infection_status'] = 'I'
        G.nodes[s]['flu_infection_step'] = 0
        G.nodes[s]['flu_ever_infected'] = True
        G.graph['flu_total_infections'] = G.graph.get('flu_total_infections', 0) + 1


def prepare_population(population_size, age_distribution, male_fraction, indigenous_fraction, mode, hiv_seeds=20, flu_seeds=20):
    """Create a population and ensure it matches the requested experimental control mode."""
    G = pop.generate_population(population_size, age_distribution, male_fraction, indigenous_fraction)

    # Your generate_population may auto-seed; normalise to the requested mode
    if mode == 'hiv_only':
        reset_flu(G)
        reset_hiv(G)
        seed_only_hiv(G, hiv_seeds)
    elif mode == 'flu_only':
        reset_hiv(G)
        reset_flu(G)
        seed_only_flu(G, flu_seeds)
    elif mode == 'both':
        # ensure both present; rebaseline so counts start from seeds below
        reset_hiv(G)
        reset_flu(G)
        seed_only_hiv(G, hiv_seeds)
        seed_only_flu(G, flu_seeds)
    else:
        raise ValueError("mode must be one of: hiv_only, flu_only, both")
    return G


# ----------------------------
# One-run simulation
# ----------------------------
def run_simulation(mode='both',
                   total_days=730,
                   relationship_formation_prob=0.1429,  # ~once/week opportunity
                   homophily=0.7,
                   min_age=16,
                   breakup_probability=0.02,
                   hiv_p_mtf=1/1234,
                   hiv_p_ftm=1/2380,
                   flu_edge_beta=0.15,
                   flu_hiv_multiplier=2.0,
                   flu_comm_contacts=5,
                   flu_comm_beta=0.1,
                   rng_seed=42):
    random.seed(rng_seed)
    np.random.seed(rng_seed)

    # ABS/Bourke defaults (from your main.py)
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
        (85, 120, 34/population_size),
    ]
    male_fraction = 1164 / population_size
    indigenous_fraction = 708 / population_size

    G = prepare_population(population_size, age_distribution, male_fraction, indigenous_fraction, mode)

    # tracking
    infected_hiv_counts = [sum(1 for _, a in G.nodes(data=True) if a.get("hiv_infection_status") == "I")]
    infected_flu_counts = [sum(1 for _, a in G.nodes(data=True) if a.get("flu_infection_status") == "I")]
    new_hiv_cases_per_day = []
    new_flu_cases_per_day = []

    # simulate
    for day in range(total_days):
        # Progress flu states regardless of mode; in hiv_only there are no flu infections so it's a no-op
        if mode in ('both', 'flu_only'):
            disease.progress_flu(G, day, incubation_period=4, infectious_period=7)

        # Relationship dynamics
        relationship.start_relationship(G, relationship_formation_prob, homophily, day, min_age, None)
        # Transmissions
        if mode in ('both', 'flu_only'):
            flu_new_today = disease.transmit_flu(
                G, day,
                edge_beta=flu_edge_beta,
                hiv_multiplier=flu_hiv_multiplier if mode == 'both' else 1.0,  # disable HIV effect in flu-only
                community_contacts=flu_comm_contacts,
                community_beta=flu_comm_beta
            )
            new_flu_cases_per_day.append(flu_new_today)

        if mode in ('both', 'hiv_only'):
            hiv_new_today = disease.transmit_hiv(
                G,
                transmission_probability_ftm=hiv_p_ftm,
                transmission_probability_mtf=hiv_p_mtf,
                current_step=day
            )
            new_hiv_cases_per_day.append(hiv_new_today)

        relationship.breakup(G, breakup_probability=breakup_probability)
        # HIV recovery / long-term status
        pop.apply_recovery(G, day)

        # Prevalence counts
        infected_hiv_counts.append(sum(1 for _, a in G.nodes(data=True) if a.get("hiv_infection_status") == "I"))
        infected_flu_counts.append(sum(1 for _, a in G.nodes(data=True) if a.get("flu_infection_status") == "I"))

    # Save a snapshot and degree
    deg = dict(G.degree())
    nx.set_node_attributes(G, deg, "degree_now")
    nx.write_graphml(G, f"bourke_{mode}_final_{total_days}d.graphml")

    # Package outputs
    outputs = {
        "G": G,
        "infected_hiv_counts": infected_hiv_counts,
        "infected_flu_counts": infected_flu_counts,
        "new_hiv_cases_per_day": new_hiv_cases_per_day,
        "new_flu_cases_per_day": new_flu_cases_per_day,
        "total_days": total_days,
        "population_size": G.number_of_nodes()
    }
    return outputs


def plot_results(mode, results):
    days = list(range(0, results["total_days"] + 1))
    popN = results["population_size"]

    # HIV prevalence (if applicable)
    if results["infected_hiv_counts"] and any(results["infected_hiv_counts"]):
        plt.figure(figsize=(8,4))
        plt.plot(days, results["infected_hiv_counts"], linestyle='-')
        plt.xlabel("Day")
        plt.ylabel("Number of HIV-infected individuals")
        plt.title(f"HIV Prevalence over Time — {mode}")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

        if results["new_hiv_cases_per_day"]:
            d = list(range(1, results["total_days"] + 1))
            cum = np.cumsum(results["new_hiv_cases_per_day"])
            per_1000 = (cum / popN) * 1000.0

            plt.figure(figsize=(9,4))
            plt.plot(d, cum, linestyle='-')
            plt.xlabel("Day")
            plt.ylabel("People ever infected")
            plt.title(f"Cumulative HIV incidence — {mode}")
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.show()

            plt.figure(figsize=(9,4))
            plt.plot(d, per_1000, linestyle='-')
            plt.xlabel("Day")
            plt.ylabel("Cumulative incidence per 1,000")
            plt.title(f"Cumulative HIV incidence per 1,000 — {mode}")
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.show()

    # Flu prevalence (if applicable)
    if results["infected_flu_counts"] and any(results["infected_flu_counts"]):
        plt.figure(figsize=(8,4))
        plt.plot(days, results["infected_flu_counts"], linestyle='-')
        plt.xlabel("Day")
        plt.ylabel("Number of Flu-infected individuals")
        plt.title(f"Influenza Prevalence over Time — {mode}")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

        if results["new_flu_cases_per_day"]:
            d = list(range(1, results["total_days"] + 1))
            cum = np.cumsum(results["new_flu_cases_per_day"])
            per_1000 = (cum / popN) * 1000.0

            plt.figure(figsize=(9,4))
            plt.plot(d, results["new_flu_cases_per_day"], linestyle='-')
            plt.xlabel("Day")
            plt.ylabel("New flu cases (incidence)")
            plt.title(f"Daily influenza incidence — {mode}")
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.show()

            plt.figure(figsize=(9,4))
            plt.plot(d, cum, linestyle='-')
            plt.xlabel("Day")
            plt.ylabel("People ever infected")
            plt.title(f"Cumulative influenza incidence — {mode}")
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.show()

            plt.figure(figsize=(9,4))
            plt.plot(d, per_1000, linestyle='-')
            plt.xlabel("Day")
            plt.ylabel("Cumulative incidence per 1,000")
            plt.title(f"Cumulative influenza incidence per 1,000 — {mode}")
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.show()


def main():
    parser = argparse.ArgumentParser(description="Run Bourke simulation in control modes: HIV-only, Flu-only, or Both.")
    parser.add_argument("--mode", choices=["hiv_only", "flu_only", "both"], default="both")
    parser.add_argument("--days", type=int, default=730)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    results = run_simulation(mode=args.mode, total_days=args.days, rng_seed=args.seed)
    print(f"Mode: {args.mode}")
    print(f"Days: {args.days}")
    print(f"Nodes: {results['population_size']}")
    print(f"HIV total infections: {results['G'].graph.get('hiv_total_infections', 0)}")
    print(f"Flu total infections: {results['G'].graph.get('flu_total_infections', 0)}")
    print("Saving GraphML to:", f"bourke_{args.mode}_final_{args.days}d.graphml")
    plot_results(args.mode, results)


if __name__ == "__main__":
    main()
