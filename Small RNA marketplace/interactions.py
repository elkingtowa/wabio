#The rule anyone samples two members of the population uniformly and independently
#nearby(pop, k), which choses one member uniformly and then chooses a second within k index elements away, to simulate interactions within a local neighborhood
#rules can be modified

def anyone(pop): 
	return random.sample(range(len(pop)), 2)

def nearby(pop, k=5): 
    i = random.randrange(len(pop))
    j = i + random.choice((1, -1)) * random.randint(1, k)
    return i, (j % len(pop))
               
def nearby1(pop): 
	return nearby(pop, 1)
