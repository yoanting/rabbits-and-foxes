
# coding: utf-8

# In[1]:

get_ipython().magic('matplotlib inline')
import numpy as np
from matplotlib import pyplot as plt


# ## Some tips and tricks for the Kinetic Monte Carlo algorithm.
# 
# For random numbers use either the `random` module or `np.random`

# In[2]:

import random
random.seed(3) # Fix seed so results don't change every time I execute
np.random.seed(3) # same thing for numpy

rates = 10*np.random.rand(4) # just to demo, use four random numbers for our event rates
rates


# In[3]:

sum_rates = sum(rates)
sum_rates


# If `sum_rates` is on average about 20 days$^{-1}$, we'd expect the average wait_time to be about 0.05 days

# In[4]:

for i in range(5):
    wait_time = random.expovariate( sum_rates )
    print(wait_time)


# To pick one of the events, it might help to generate a uniform number between 0 and `sum_rates`

# In[5]:

for i in range(5):
    choice = random.uniform(0, sum_rates)
    print(choice)


# ## Tricks for slicing arrays

# To demonstrate this we need to generate arrays of foxes and times, so we'll quickly do the ODE solve version:

# In[6]:

from scipy.integrate import odeint
k1 = 0.015
k2 = 0.00004
k3 = 0.0004
k4 = 0.04
end_time = 600.
def rates(variables, time):
    """
    Return the right hand side of the ODE
    """
    rabbits, foxes = variables
    rate_rabbits = (k1 * rabbits - k2 * rabbits * foxes)
    rate_foxes = (k3 * rabbits * foxes - k4 * foxes)
    return (rate_rabbits, rate_foxes)

times = np.arange(0, end_time)
initial_conditions = (400., 200.)
result = odeint(rates, initial_conditions, times)
rabbits = result[:,0]
foxes = result[:,1]
plt.plot(times, rabbits, label='rabbits')
plt.plot(times, foxes, label='foxes')
plt.legend(loc="best") # put the legend at the best location to avoid overlapping things
plt.show()
print("Second peak occurs at",(foxes.argmax(),foxes.max()))


# First we will use slicing to select every 10th element from the times, rabbits, and foxes arrays
# and store in arrays called `times_`, `rabbits_`, `foxes_` just so that we can demonstrate a few tricks without printing 600 numbers that are hard to read.

# In[7]:

times_, rabbits_, foxes_  = times[::10], rabbits[::10], foxes[::10]
times_, rabbits_, foxes_


# If we call a comparison operator on an array, we get an array of Boolean values like this

# In[8]:

times_ > 200


# If we multiply a number by a boolean, then the True is treated as 1 and the False as 0, like this:

# In[9]:

5 * (times_>200)


# We can chain these together, for example to make an array of the number of foxes but only containing values corresponding to times > 200 and foxes > 100:

# In[10]:

foxes_ * (times_>200) * (foxes_>100)


# We just sliced out where we want to look for a second peak:

# In[11]:

plt.plot(times,foxes*(times>200)*(foxes>100))
plt.show()


# Now if we want to find the index of the highest of those values, we can use `np.argmax()` (check its documentation if unsure)

# In[12]:

i = np.argmax(foxes_*(times_>200)*(foxes_>100))
i


# then to find the location (on a fox-v-time plot) of that peak

# In[13]:

times_[i], foxes_[i]


# These tricks relied upon using numpy arrays rather than lists.
# If you know how big something is going to be before you start, then you can create an empty numpy array to store it in, eg. using `np.zeros()`

# In[14]:

runs = 10
average_peak_time = np.zeros(runs)
for i in range(runs):
    average_peak_time[i] = random.random()
average_peak_time


# But if you're not sure how long something will be then start with an empty list, append to it, then convert it into a numpy array when you're finished:

# In[15]:

times = []
foxes = []
time = 0
while time < 600:
    time += wait_time
    times.append(time)
    foxes.append(1)
times = np.array(times)
foxes = np.array(foxes)
times, foxes

