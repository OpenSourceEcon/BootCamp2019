import gym
from gym import spaces
from gym.utils import seeding
import numpy as np
import math

alpha_param = 0.3
delta_param = 0.1
klow_param = 0.1
khigh_param = 4
num_action_param = 15
num_ex_shocks = 2
trans_matrix = np.array([[0.75, 0.25], [0.15, 0.85]])
tfp = np.array([0.9, 1.1])




class Ramsey(gym.Env):
    def __init__(self):
        self.alpha = alpha_param
        self.delta = delta_param
        
        self.klow = klow_param

        self.khigh = khigh_param
        
        self.num_actions = num_action_param
        
        self.klist=np.linspace(self.klow, self.khigh, self.num_actions)
        
        self.action_space = spaces.Discrete(self.num_actions)
        self.observation_space = spaces.MultiDiscrete([num_ex_shocks, num_action_param])
        
        self.observation = [0, 0]
        
        self.seed()
        self.reset()
        
    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]
    
    def set_observation(self, obs):
        assert self.observation_space.contains(obs), 'observation not in observation set'
        
        self.observation = obs
        
    def get_klist(self):
        return self.klist
    
    def step(self, action):
        assert self.action_space.contains(action), 'action not in observation set'
        cur_shock = self.observation[0]
        probs_next = trans_matrix[cur_shock, :]

        cur_k = self.klist[self.observation[1]]
        next_k = self.klist[action]
        
        cur_cons = tfp[cur_shock] * cur_k ** self.alpha + (1 - self.delta) * cur_k - next_k
        
        cons_thresh = 1E-6

        if cur_cons > cons_thresh:

            reward = np.log(cur_cons)
            done =  False
        
        else:
            reward = np.log(cons_thresh) + (1. / cons_thresh) * (cur_cons - cons_thresh)

            done = True
        
        rn = np.random.rand()
        for i in range(num_ex_shocks):
            if rn <= np.sum(probs_next[0:i+1]):
                next_shock = i
                break


        self.observation = np.array([next_shock, action], dtype = np.int8)
        
        return self.observation, reward, done
    
    def get_legal_actions(self):
        
        legal_actions = []
        #print('legal_actions', legal_actions)
        cur_k = self.klist[self.observation]
        
        
        for i in range(self.num_actions):
            next_k = self.klist[i]
            cur_cons = cur_k ** self.alpha + (1 - self.delta) * cur_k - next_k
            if cur_cons > 0:
                legal_actions.append(i)
            else:
                return legal_actions
        return legal_actions
    
            
    
    def reset(self, start_obs = 0):
        self.observation = start_obs
        return self.observation
