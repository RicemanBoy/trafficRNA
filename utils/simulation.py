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
        self.tf = 0                        #0 = no dynamic k_on, 1 = dynamic k_on

        self.k_on = tf_probs            #original value of the k_on rates
        self.k_off = tf_probs

        self.box_p = np.ones((length,))*tf_probs        #probabilities for tf to dock (and undock)

        self.meth_p = m_on        #probability for methylization of the boxes

        self.meth_off = m_off           #probability for methylazation to go away

        self.meth_p_orig = m_on         #backup of the original_values

        self.production = []            #tracks how much protein is produced each step

        self.track = [0 for i in range(length)]    #tracks which boxes are occupied

        # #TEST to reach equil faster##########
        # for i, val in enumerate(self.track):
        #     r, r1, r2, r3, r4 = np.random.rand(), np.random.rand(), np.random.rand(), np.random.rand(), np.random.rand()
        #     if r < 0.5:
        #         if r1<tf_probs:
        #             self.track[i] = self.boxes[i]
        #         elif r2<self.meth_p:
        #             self.track[i] = "methyl"
        #     else:
        #         if r3<self.meth_p:
        #             self.track[i] = "methyl"
        #         elif r4<tf_probs:
        #             self.track[i] = self.boxes[i]
        # self.dyn_methy()
        # ####################################

    def count_TF(self):                #counts how many boxes are filled with TFs
        occupied = sum(1 for i in self.track if i != "methyl" and i != 0)
        return occupied
    
    def protein(self):              #returns the amount of protein in this timestep
        filled = self.count_TF()
        return filled

    def dyn_methy(self):                #change m_on based on the total methylation
        if self.methyl == 1:
            count = self.track.count("methyl")
            self.meth_p = self.meth_p_orig*(1+(1/self.L)*count)
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

            elif self.track[pos] == self.boxes[pos] and random3 <= self.k_off:          #TF --> 0
                self.track[pos] = 0
            
            elif self.track[pos] == "methyl" and random4 <= self.meth_off:         #Methyl --> 0
                self.track[pos] = 0
        else:
            if self.track[pos] == 0 and random2 <= self.meth_p:          #0 --> Methyl
                self.track[pos] = "methyl"
            
            elif self.track[pos] == 0 and random1 <= self.box_p[pos]:          #0 --> TF
                self.track[pos] = self.boxes[pos]

            elif self.track[pos] == self.boxes[pos] and random3 <= self.k_off:          #TF --> 0
                self.track[pos] = 0
            
            elif self.track[pos] == "methyl" and random4 <= self.meth_off:         #Methyl --> 0
                self.track[pos] = 0
        
        self.dyn_methy()        #update dynamic methylation

class simulation:
    def __init__(self, timesteps: int):
        self.timesteps = timesteps
        self.timeline = [i for i in range(timesteps)]

    def dyn_tf(self, gene1: gene, gene2: gene):
        if gene1.tf == 1 and gene2.tf == 1:
            for i in range(gene1.L):
                if gene2.track[i] == letters[i]:
                    gene1.box_p[i] = gene1.k_on*0.5
                else:
                    gene1.box_p[i] = gene1.k_on

                if gene1.track[i] == letters[i]:
                    gene2.box_p[i] = gene2.k_on*0.5
                else:
                    gene2.box_p[i] = gene2.k_on

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
            self.dyn_tf(gene1=gene1, gene2=gene2)        #dynamic rates for TF, boxes are indepent of each other for now
            gene1.time += 1
            gene2.time += 1
            
        
    def average(self, length: int, tf_probs: float, m_on: float, m_off: float, shots: int, m_switch = True, tf_switch = True):              #run the simulation #shots times and return the avg/std
        avg1, avg2, std1, std2 = [],[],[],[]
        corr = []

        rate = min(tf_probs, m_on, m_off)
        equil = int(np.log(0.01)/np.log(1-rate))   #determine equil based on smallest rate, set the threshhold to reach equil 0.01 as fixed, can be changed

        timewindow, tw1, tw2, tw3 = 0.5, 1, 2, 4            #with respect to equil

        corrx = np.array([0 for i in range(int(self.timesteps-(1+timewindow)*equil))])
        corrx1 = np.array([0 for i in range(int(self.timesteps-(1+tw1)*equil))])
        corrx2 = np.array([0 for i in range(int(self.timesteps-(1+tw2)*equil))])
        corrx3 = np.array([0 for i in range(int(self.timesteps-(1+tw3)*equil))])

        total1, total2 = np.array([0 for i in range(self.timesteps)]), np.array([0 for i in range(self.timesteps)])

        for z in range(shots):
            gene1, gene2 = gene(length, tf_probs, m_on, m_off), gene(length, tf_probs, m_on, m_off)
            if m_switch:                    #decide if neighbours affect methlylzation
                gene1.methyl = 1
                gene2.methyl = 1
            if tf_switch:
                gene1.tf = 1
                gene2.tf = 1
            self.run(gene1, gene2)

            total1 = np.vstack((total1, gene1.production))
            total2 = np.vstack((total2, gene2.production))
      
            try:
                if np.isnan(np.corrcoef(gene1.production[equil:], gene2.production[equil:])[0][1]):
                    corr.append(0)
                else:
                    corr.append(np.corrcoef(gene1.production[equil:], gene2.production[equil:])[0][1])
            except:
                try:
                    equil_ = int(equil/2)
                    if np.isnan(np.corrcoef(gene1.production[equil_:], gene2.production[equil_:])[0][1]):
                        corr.append(0)
                    else:
                        corr.append(np.corrcoef(gene1.production[equil_:], gene2.production[equil_:])[0][1])
                except:
                    equil_ = int(equil_/2)
                    if np.isnan(np.corrcoef(gene1.production[equil_:], gene2.production[equil_:])[0][1]):
                        corr.append(0)
                    else:
                        corr.append(np.corrcoef(gene1.production[equil_:], gene2.production[equil_:])[0][1])
            
            corr_idk = np.array([0])
            for i in range(int(self.timesteps-(1+timewindow)*equil)):
                corr_piece = np.corrcoef(gene1.production[equil+i:int((1+timewindow)*equil)+i], gene2.production[equil+i:int((1+timewindow)*equil)+i])[0][1]
                if np.isnan(corr_piece):
                    corr_idk = np.append(corr_idk, 0)
                else:
                    corr_idk = np.append(corr_idk, corr_piece)
            corr_idk = np.delete(corr_idk, 0)
            corrx = np.vstack((corrx, corr_idk))

            corr_idk = np.array([0])
            for i in range(int(self.timesteps-(1+tw1)*equil)):
                corr_piece = np.corrcoef(gene1.production[equil+i:int((1+tw1)*equil)+i], gene2.production[equil+i:int((1+tw1)*equil)+i])[0][1]
                if np.isnan(corr_piece):
                    corr_idk = np.append(corr_idk, 0)
                else:
                    corr_idk = np.append(corr_idk, corr_piece)
            corr_idk = np.delete(corr_idk, 0)
            corrx1 = np.vstack((corrx1, corr_idk))

            corr_idk = np.array([0])
            for i in range(int(self.timesteps-(1+tw2)*equil)):
                corr_piece = np.corrcoef(gene1.production[equil+i:int((1+tw2)*equil)+i], gene2.production[equil+i:int((1+tw2)*equil)+i])[0][1]
                if np.isnan(corr_piece):
                    corr_idk = np.append(corr_idk, 0)
                else:
                    corr_idk = np.append(corr_idk, corr_piece)
            corr_idk = np.delete(corr_idk, 0)
            corrx2 = np.vstack((corrx2, corr_idk))

            corr_idk = np.array([0])
            for i in range(int(self.timesteps-(1+tw3)*equil)):
                corr_piece = np.corrcoef(gene1.production[equil+i:int((1+tw3)*equil)+i], gene2.production[equil+i:int((1+tw3)*equil)+i])[0][1]
                if np.isnan(corr_piece):
                    corr_idk = np.append(corr_idk, 0)
                else:
                    corr_idk = np.append(corr_idk, corr_piece)
            corr_idk = np.delete(corr_idk, 0)
            corrx3 = np.vstack((corrx3, corr_idk))

        corrx, corrx1, corrx2, corrx3 = np.delete(corrx, (0), axis=0), np.delete(corrx1, (0), axis=0), np.delete(corrx2, (0), axis=0), np.delete(corrx3, (0), axis=0)
        total1 = np.delete(total1, (0), axis=0)
        total2 = np.delete(total2, (0), axis=0)
        for i in total1.T:
            avg1.append(np.mean(i))
            std1.append(np.std(i))
        for i in total2.T:
            avg2.append(np.mean(i))
            std2.append(np.std(i))

        return avg1, std1, avg2, std2, corr, corrx, corrx1, corrx2, corrx3

    




