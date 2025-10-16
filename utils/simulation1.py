import numpy as np

class rna:
    def __init__(self, L, tau, steps):
        self.L = L
        self.tau = tau
        self.steps = steps

        self.track = [0 for i in range(L)]      #tracks, where the ribosomes are

        self.time = 0           #for every timestep --> +tau

    def time_update(self):
        self.time += self.tau
    
    def occupation(self):           #returns where the ribosomes are
         return self.track

class ribo:
    def __init__(self, rna):        #whenever creating a ribo, check for rna if spot is empty!
            self.pos = 0
            rna.track[0] += 1
    
    def move(self, rna):
        if self.pos+1 == rna.L:
             self.pos = "protein"
             rna.track[self.pos] -= 1

        if rna.track[self.pos+1] == 0:
            rna.track[self.pos] -= 1
            rna.track[self.pos+1] += 1
            self.pos += 1

    def position(self):
        return self.pos

