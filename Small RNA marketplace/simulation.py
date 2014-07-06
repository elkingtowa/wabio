#simulation of an economic marketplace in which there is a population of actors, each of which has a level of 
#wealth (a single number) that changes over time. On each time step two agents (chosen by an interaction rule) 
#interact with each other and exchange wealth (according to a transaction rule)

#The function simulate does the work; it runs the interaction function to find two actors, 
#then calls the transaction function to figure out how to split their wealth, and repeats this T times. 
#The only other thing it does is record results. Every so-many steps, it records some summary statistics of 
#the population (by default, this will be every 25 steps)

#small rna marketplace

def simulate(population, transaction_fn, interaction_fn, T, percentiles, record_every):
    "Run simulation for T steps; collect percentiles every 'record_every' time steps."
    results = []
    for t in range(T):
        i, j = interaction_fn(population)
        population[i], population[j] = transaction_fn(population[i], population[j]) 
        if t % record_every == 0:
            results.append(record_percentiles(population, percentiles))
    return results

def report(distribution=gauss, transaction_fn=random_split, interaction_fn=anyone, N=N, mu=mu, T=5*N, 
           percentiles=(1, 10, 25, 33.3, 50, -33.3, -25, -10, -1), record_every=25):
    "Print and plot the results of the simulation running T steps." 
    # Run simulation
    population = sample(distribution, N, mu)
    results = simulate(population, transaction_fn, interaction_fn, T, percentiles, record_every)
    # Print summary
    print('Simulation: {} * {}(mu={}) for T={} steps with {} doing {}:\n'.format(
          N, name(distribution), mu, T, name(interaction_fn), name(transaction_fn)))
    fmt = '{:6}' + '{:10.2f} ' * len(percentiles)
    print(('{:6}' + '{:>10} ' * len(percentiles)).format('', *map(percentile_name, percentiles)))
    for (label, nums) in [('start', results[0]), ('mid', results[len(results)//2]), ('final', results[-1])]:
        print fmt.format(label, *nums)
    # Plot results
    for line in zip(*results):
        plt.plot(line)
    plt.show()

def record_percentiles(population, percentiles):
    "Pick out the percentiles from population."
    population = sorted(population, reverse=True)
    N = len(population)
    return [population[int(p*N/100.)] for p in percentiles]

def percentile_name(p):
    return ('median' if p == 50 else 
            '{} {}%'.format(('top' if p > 0 else 'bot'), abs(p)))
    
def name(obj):
    return getattr(obj, '__name__', str(obj))