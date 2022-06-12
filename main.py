import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class Configuracoes:
    def __init__(self):
        self.width = 40     # cm
        self.height = 80     # cm

        self.fc = 280
        self.fy = 4200
        self.ec_max = 0.003                     # concrete maximum strain (Constant)
        self.es_min = 0.005                     # 受拉鋼筋之淨拉應變
        self.Ec = 15000 * self.fc ** 0.5
        self.Es = 2.04 * 10 ** 6                # steel modulus (kgf/cm2)
        self.ey = self.fy / self.Es             # steel strain at yield point

        self.rebars = [{"area":12,"coord_x":20,"coord_y":7.2},
                       {"area":12,"coord_x":20,"coord_y":72.8}]

        self.rebars_lowest = min(rebar["coord_y"] for rebar in self.rebars)
        self.rebars_totalArea = sum(rebar['area'] for rebar in self.rebars)


        self.beta1 = self.__beta1()


    def __beta1(self):
        '''
            # Beta1 (whitney stress block)
        '''
        if self.fc <= 280:
            beta1 = 0.85
        else:
            beta1 = max(0.85 - 0.05 * ((self.fc - 280) / 70), 0.65)
        return beta1

    def steel_strain_stress_at_depth(self, neutral_axis_depth, depth_of_interest):
        # 拉力為正；壓力為負
        strain = -self.ec_max/neutral_axis_depth * (neutral_axis_depth - depth_of_interest)
        if strain > 0:
            stress = min(strain * self.Es, +self.fy)
        else:
            stress = max(strain * self.Es, -self.fy)
        return strain, stress

    def strength_factor(self, epsilon_t, is_spiral = False):
        '''
        # Strength reduction factor
        '''
        # 橫箍筋 0.65 ；螺箍筋 0.7 
        phi = 0.65 if is_spiral == False else 0.7
        epsilon_y = self.fy / self.Es
        if epsilon_t >= 0.005:
            phi, classify = 0.9, "Tension-controlled"
        elif epsilon_t <= epsilon_y:
            phi, classify = phi, "Compression-controlled"
        else:
            phi_factor = phi + 0.25 * (epsilon_t - self.ey) / (0.005 - self.ey)
            phi, classify = phi_factor, "Transition"    
        return phi, classify

    def forces_moments(self, neutral_axis_depth):
        
        a = self.beta1 * neutral_axis_depth

        # area_sic : area of steel rebar in concrete compress area
        steel_area_in_compress_area = 0
        steel_axial_forces = 0
        steel_moments = 0
        for rebar in self.rebars :
            # multiplying the stress in each layer by the area in each layer.
            strain, stress = self.steel_strain_stress_at_depth(neutral_axis_depth, rebar["coord_y"])
            steel_axial_forces += stress * rebar['area']
            steel_moments += stress * rebar['area'] * rebar["coord_y"]
            
            # 混凝土壓力區內之鋼筋面積        
            if rebar["coord_y"] < a:
                steel_area_in_compress_area += rebar['area']
        
        concrete_axial_forces = 0.85 * self.fc * (a * self.width - steel_area_in_compress_area)
        Pn = concrete_axial_forces - steel_axial_forces
        Mn = -concrete_axial_forces * a / 2 + steel_moments + Pn*self.height/2
        Pu = phi * Pn
        Mu = phi * Mn
        ecc = Mn / Pn  
        return ecc, Pn, Pu, Mn, Mu

class Rebars:
    def __init__(self, rebars_obj):
        self.rebars = rebars_obj
        self.lowest = self.__lowest()
        self.total_area = self.__total_area()
        self.coords = self.__coords()

    def __lowest(self):
        coord_ymin = min(rebar["coord"][1] for rebar in self.rebars)
        return coord_ymin 

    def __total_area(self):
        area = sum(rebar['area'] for rebar in self.rebars)
        return area

    def __coords(self):
        return [rebar["coord"] for rebar in self.rebars]   

class Rec:
    pass

# 鋼筋資料
rebars = [{"area":6,"coord":( 5,  5)},
          {"area":6,"coord":( 5, 75)},
          {"area":6,"coord":(35,  5)},
          {"area":6,"coord":(35, 75)}]
rb = Rebars(rebars)

cfg = Configuracoes()
nom_load = []
ult_load = []
nom_moment = []
ult_moment = []
eccentricity = []
phi_factor = []


c_value = np.arange(cfg.height,8,-0.5)
for c in c_value:
    es, fs = cfg.steel_strain_stress_at_depth(c, cfg.height-cfg.rebars_lowest )
    phi, classify =  cfg.strength_factor(es)
    ecc, Pn, Pu, Mn, Mu = cfg.forces_moments(c)

    if ecc <= 1.5 * cfg.height:
        nom_load.append(round(Pn))
        ult_load.append(round(Pu))
        nom_moment.append(round(Mn))
        ult_moment.append(round(Mu))
        eccentricity.append(round(ecc))
        phi_factor.append(round(phi,2))

dict = {"ecc": eccentricity,
        "Pn": nom_load,
        "Pu": ult_load,
        "Mn": nom_moment,
        "Mu": ult_moment}

df = pd.DataFrame(dict)

print(df)

sns.set_style("darkgrid")
fig, ax = plt.subplots(figsize=(13,7))
ax = sns.scatterplot(x= "Mn", y= "Pn" , data =df, color="g", label="Normal")
sns.scatterplot(x= "Mu", y= "Pu" , data =df, color="r", label="Ultimate")
ax.set_xlabel("Moment Mn")
ax.set_ylabel("Axial Load Pn")
plt.title("P-M Interation Diagram")
plt.show()

print(rb.total_area)
print(rb.lowest)
print(rb.coords)