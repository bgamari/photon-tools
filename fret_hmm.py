from collections import namedtuple
import random
from numpy import array
from numpy.random import random_sample, randint, poisson
#from hmm import HMM, baum_welsh, viterebi

def weighted_choice(choices, probs):
        #assert sum(probs) == 1
        r = random.random()
        for c,p in zip(choices, probs):
                if r < p: return c
                r -= p
        raise RuntimeError()

Model = namedtuple('Model', 'n_states n_obs start_prob trans_prob emission')
def random_model(n_states, n_obs):
        start_prob = random_sample(n_states)
        start_prob /= sum(start_prob)

        emission = randint(0, 150, (n_states, n_obs))

        transition_prob = []
        for i in xrange(n_states):
                stay_prob = random.random()
                t = stay_prob * random_sample(n_states)
                t[i] = stay_prob
                t /= sum(t)
                transition_prob.append(t)

        transition_prob = array(transition_prob)
        return Model(n_states, n_obs, start_prob, transition_prob, emission)

def random_data(model, length, noise=True):
        state_seq = []
        state = weighted_choice(xrange(model.n_states), model.start_prob)
        data = []
        for i in xrange(length):
                state_seq.append(state)
                datum = model.emission[state]
                if noise:
                        datum = poisson(datum) 
                data.append(datum)
                state = weighted_choice(xrange(model.n_states), model.trans_prob[state])

        return data, state_seq


n_states = 6
model = random_model(n_states, 1)
data, seq = random_data(model, 100)

from matplotlib import pyplot as pl
pl.plot(data)
pl.plot(seq)
pl.show()

hmm = HMM(n_states)
baum_welsh(hmm, data, graph=True)

