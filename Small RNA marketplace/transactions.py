#ransaction rules will decide how to split up the pot of X+Y total wealth

def random_split(X, Y):
    "Take all the money in the pot and divide it randomly between X and Y."
    pot = X + Y
    m = random.uniform(0, pot)
    return m, pot - m
    
def winner_take_most(X, Y, most=3/4.): 
    "Give most of the money in the pot to one of the parties."
    pot = X + Y
    m = random.choice((most * pot, (1 - most) * pot))
    return m, pot - m

def winner_take_all(X, Y): 
    "Give all the money in the pot to one of the actors."
    return winner_take_most(X, Y, 1.0)

def redistribute(X, Y): 
    "Give 55% of the pot to the winner; 45% to the loser."
    return winner_take_most(X, Y, 0.55)

def split_half_min(X, Y):
    """The poorer actor only wants to risk half his wealth; 
    the other actor matches this; then we randomly split the pot."""
    pot = min(X, Y)
    m = random.uniform(0, pot)
    return X - pot/2. + m, Y + pot/2. - m