#!/usr/bin/python

from collections import namedtuple
import random
from numpy import array, sum, mean, std
from numpy.random import random_sample, randint, poisson
import ghmm
from matplotlib import pyplot as pl

Model = namedtuple('Model', 'n_states n_obs start_prob trans_prob emissions')

def weighted_choice(choices, probs):
        #assert sum(probs) == 1
        r = random.random()
        for c,p in zip(choices, probs):
                if r < p: return c
                r -= p
        raise RuntimeError()

def random_model(n_states, n_obs):
        """ Create a ranodm model with the given number of states and observables """
        emissions = randint(0, 1500, (n_states, n_obs))
        return random_model_from_emissions(emissions)

def random_model_from_emissions(emissions):
        """ Create a random model using a given set of emission parameters.
            emissions should be an MxN matrix (M=number of states, N=number of observables) """
        n_states, n_obs = emissions.shape
        start_prob = random_sample(n_states)
        start_prob /= sum(start_prob)

        transition_prob = []
        for i in xrange(n_states):
                stay_prob = random.random()
                t = stay_prob * random_sample(n_states)
                t[i] = stay_prob
                t /= sum(t)
                transition_prob.append(t)

        transition_prob = array(transition_prob)
        return Model(n_states, n_obs, start_prob, transition_prob, emissions)

def random_data(model, length, noise=True):
        """ Generate emissions data from the given model """
        state_seq = []
        state = weighted_choice(xrange(model.n_states), model.start_prob)
        data = []
        for i in xrange(length):
                state_seq.append(state)
                datum = model.emissions[state]
                if noise:
                        datum = poisson(datum) 
                data.append(datum)
                state = weighted_choice(xrange(model.n_states), model.trans_prob[state])

        return data, state_seq

def transition_matrix(hmm):
        """ Get the transition matrix from a ghmm.HMM as a numpy array """
        get_row = lambda row : [ hmm.getTransition(row,col) for col in range(hmm.N) ]
        return array([ get_row(row) for row in range(hmm.N) ])

# Generate a model to pull data from
n_states = 6
model = random_model(n_states, 1)
print "Model:"
print model.trans_prob[0,:]

# Optional: Plot a sample of data
if False:
        data, seq = random_data(model, 100000)
        pl.plot(data)
        pl.plot(seq)
        pl.show()
        pl.clf()

# Optional: Plot FPT distribution
if True:
        data, seq = random_data(model, 100000)
        pl.suptitle("Original Model")
        pl.hist(data, 100)
        text = "State Emissions:\n" + '\n'.join( ['%d:  %d' % (s, model.emissions[s]) for s in range(n_states)] )
        pl.figtext(0.7, 0.6, text)
        pl.savefig('model.png')
        pl.clf()

# Generate a data set with which to track our convergence
# This will not be learned from, only tested for likelihood
dom = ghmm.Float()
data, seq = random_data(model, 100000)
test_data = ghmm.EmissionSequence(dom, data)

# Try teaching several randomly initialized models
print
print "Learn:"
for i in range(5):
        # Setup HMM with a new random model
        new = random_model(n_states, 1)
        B = [ [float(e[0]), float(e[0])] for e in model.emissions ]  # mu, sigma
        hmm = ghmm.HMMFromMatrices(dom, ghmm.GaussianDistribution(dom),
                                   new.trans_prob, B, new.start_prob)

        if False:
                data = hmm.sampleSingle(100000)
                pl.suptitle("Initial Training Model %d" % i)
                pl.hist(data, 100)
                pl.savefig('initial-%d.png' % i)
                pl.clf()

        tm = transition_matrix(hmm)
        print tm[0,:], '%e' % mean((tm - model.trans_prob)**2), '%e' % hmm.loglikelihoods(test_data)[0]

        # Training iterations
        for j in range(10):
                data, seq = random_data(model, 100000)
                seq = ghmm.EmissionSequence(dom, data)
                hmm.baumWelch(seq, 50, 0.1)

                if False:
                        tm = transition_matrix(hmm)
                        print tm[0,:], '%e' % mean((tm - model.trans_prob)**2), '%e' % hmm.loglikelihoods(test_data)[0]

        tm = transition_matrix(hmm)
        print tm[0,:], '%e' % mean((tm - model.trans_prob)**2), '%e' % hmm.loglikelihoods(test_data)[0]
        print

        if True:
                data = hmm.sampleSingle(100000)
                pl.suptitle("Training Model %d" % i)
                pl.hist(data, 100)
                pl.savefig('trained-%d.png' % i)
                pl.clf()

pl.show()
