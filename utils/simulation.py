#new file based on the new concept, what Lizanza and me discussed

import numpy as np
import regex

letters = ''.join(regex.findall(r'\p{L}', ''.join(chr(i) for i in range(0x110000))))

class gene:
    def __init__(self, length: int, tf_probs: float, m_on: float, m_off: float):
        assert length in range(141028), "Exception message raised because parameter 'length' is larger than 51!"
        self.boxes = [letters[i] for i in range(length)]        #define the boxes for the transcription factors
        self.L = length                                                 #each box has its own TF

        self.time = 0

        self.methyl = 0                 #0 = neighbours dont affect p_m, 1 = neighbours DO affect p_m

        #self.box_p = np.random.rand(length)*tf_probs        #probabilities for tf to dock (and undock)
        #self.meth_p = np.random.rand(length)*m_probs        #probability for methylization of the boxes (and also unmethylization)

        self.box_p = np.ones((length,))*tf_probs        #probabilities for tf to dock (and undock)
        self.meth_p = m_on        #probability for methylization of the boxes (and also unmethylization)

        self.meth_off = m_off           #probability for methylazation to go away

        self.meth_p_orig = m_on         #backup of the original_values

        self.production = []            #tracks how much protein is produced each step

        self.track = [0 for i in range(length)]    #tracks which boxes are occupied

    def count_TF(self):                #counts how many boxes are filled with TFs
        occupied = sum(1 for i in self.track if i != "methyl" and i != 0)
        return occupied
    
    def protein(self):              #returns the amount of protein in this timestep
        filled = self.count_TF()
        return filled

    def dyn_methy(self):                #change m_on based on the total methylation
        if self.methyl == 1:
            count = self.track.count("methyl")
            self.meth_p = self.meth_p_orig*(1+0.2*count)
        # if self.L > 1:
        #     if pos == 0:
        #         if self.track[1] == "methyl":
        #             count += 1
        #     elif pos == self.L - 1:
        #         if self.track[self.L-2] == "methyl":
        #             count += 1
        #     else:
        #         if self.track[pos-1] == "methyl":
        #             count += 1
        #         if self.track[pos+1] == "methyl":
        #             count += 1
        #     if self.methyl == 1:
        #         self.meth_p[pos] = self.meth_p_orig[pos]*(1+count)

    def methyl_turn(self, click: bool):         #turn on/off dyn. methylation
        if click:
            self.methyl = 1
        else:
            self.methyl = 0


    def update(self, pos: int):

        random1, random2, random3, random4 = np.random.rand(), np.random.rand(), np.random.rand(), np.random.rand()
        rand_upd = np.random.rand()

        if rand_upd <= 0.5:     #Randomize if first check for TF OR Methyl
            if self.track[pos] == 0 and random1 <= self.box_p[pos]:          #0 --> TF
                self.track[pos] = self.boxes[pos]
            
            elif self.track[pos] == 0 and random2 <= self.meth_p:          #0 --> Methyl
                self.track[pos] = "methyl"

            elif self.track[pos] == self.boxes[pos] and random3 <= self.box_p[pos]:          #TF --> 0
                self.track[pos] = 0
            
            elif self.track[pos] == "methyl" and random4 <= self.meth_off:         #Methyl --> 0
                self.track[pos] = 0
        else:
            if self.track[pos] == 0 and random2 <= self.meth_p:          #0 --> Methyl
                self.track[pos] = "methyl"
            
            elif self.track[pos] == 0 and random1 <= self.box_p[pos]:          #0 --> TF
                self.track[pos] = self.boxes[pos]

            elif self.track[pos] == self.boxes[pos] and random3 <= self.box_p[pos]:          #TF --> 0
                self.track[pos] = 0
            
            elif self.track[pos] == "methyl" and random4 <= self.meth_off:         #Methyl --> 0
                self.track[pos] = 0
        
        self.dyn_methy()        #update dynamic methylation

class simulation:
    def __init__(self, timesteps: int):
        self.timesteps = timesteps
        self.timeline = [i for i in range(timesteps)]

    def run(self, gene1: gene, gene2: gene):
        all_boxes = gene1.L + gene2.L
        order = np.random.permutation(all_boxes)
        len = gene1.L

        for k in range(self.timesteps):
            gene1.production.append(gene1.protein())
            gene2.production.append(gene2.protein())
            for i in order:
                if i - len < 0:
                    gene1.update(i)
                else:
                    gene2.update(i-len)
            gene1.time += 1
            gene2.time += 1
            
        
    def average(self, length: int, tf_probs: float, m_on: float, m_off: float, shots: int, m_switch = True):              #run the simulation #shots times and return the avg/std
        avg1, avg2, std1, std2 = [],[],[],[]
        corr = []
        total1, total2 = np.array([0 for i in range(self.timesteps)]), np.array([0 for i in range(self.timesteps)])

        for z in range(shots):
            gene1, gene2 = gene(length, tf_probs, m_on, m_off), gene(length, tf_probs, m_on, m_off)
            if m_switch:                    #decide if neighbours affect methlylzation
                gene1.methyl = 1
                gene2.methyl = 1
            self.run(gene1, gene2)

            total1 = np.vstack((total1, gene1.production))
            total2 = np.vstack((total2, gene2.production))

            corr.append(np.corrcoef(gene1.production[200:], gene2.production[200:])[0][1])

        total1 = np.delete(total1, (0), axis=0)
        total2 = np.delete(total2, (0), axis=0)
        for i in total1.T:
            avg1.append(np.mean(i))
            std1.append(np.std(i))
        for i in total2.T:
            avg2.append(np.mean(i))
            std2.append(np.std(i))

        return avg1, std1, avg2, std2, corr

    




