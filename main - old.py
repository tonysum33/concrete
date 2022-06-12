import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class Configuracoes:
    def __init__(self):
        self.width = 40     # cm
        self.height = 80     # cm
        self.cover = 7.2
        self.deff = self.height - self.cover
        self.dprime = 7.2

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

        self.Ast = 12
        self.Asc = 12

        self.beta1 = self.__beta1()

        # self.max_rebar_tension_strain = 0.004
        # self.phi = 0.9
        # self.is_seismic_zone = True


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
        T = self.Ast * min(self.fy,fs)
        Cc = 0.85 * self.fc * a * self.width
        Cs = self.Asc * min(self.fy, fsprime) - 0.85 * self.fc
        Pn = Cc + Cs - T
        Mn = (Cc *(self.deff - a/2) + Cs * (self.deff -self.dprime))   - (self.height/2 - self.dprime) * Pn
        Pu = phi * Pn
        Mu = phi * Mn
        ecc = Mn / Pn  
        return ecc, Pn, Pu, Mn, Mu

class Rebars:
    def __init__(self):
        self.rebars = [{"area":6,"coord":( 5,  5)},
                       {"area":6,"coord":( 5, 75)},
                       {"area":6,"coord":(35,  5)},
                       {"area":6,"coord":(35, 75)}]

        self.lowest = self.__lowest()
        self.total_area = self.__total_area()

    def __lowest(self):
        coord_ymin = min(rebar["coord"][1] for rebar in self.rebars)
        return coord_ymin 

    def __total_area(self):
        area = sum(rebar['area'] for rebar in self.rebars)
        return area

# class Rcdesign:
#     def __init__(self, cfg: Configuracoes):
#         """
#         矩形梁鋼筋設計
#         :type mat: object
#         """

#         self.fc = cfg.fc
#         self.fy = cfg.fy
#         self.Es = cfg.Es
#         self.beta1 = cfg.beta1
#         self.bw = cfg.width
#         self.depth = cfg.height
#         self.phi = cfg.phi

#         self.cover = cfg.cover
#         self.epsilon_c_max = cfg.ec_max
#         self.epsilon_s_min = cfg.es_min
#         self.is_seismic_zone = cfg.is_seismic_zone
#         self.deff = cfg.deff
#         self.ps_min = self.__ps_min()
#         self.ps_max = self.__ps_max()
#         self.as_min = self.ps_min * self.bw * self.deff
#         self.as_max = self.ps_max * self.bw * self.deff

#     def moment_to_rebar(self, mu: float):
#         """
#         :param mu: 係數化彎矩 kgf-cm
#         :return:
#         """
#         # a : the depth of the compression block
#         a = self.deff - (self.deff ** 2 - 2 * abs(mu) / (0.85 * self.fc * self.phi * self.bw)) ** 0.5
#         c_max = self.epsilon_c_max / (self.epsilon_c_max + self.epsilon_s_min) * self.deff
#         a_max = self.beta1 * c_max

#         if a <= a_max:
#             # 單筋梁
#             as_t_req = mu / (self.phi * self.fy * (self.deff - a / 2))
#             as_c_req = 0

#         else:
#             # 雙筋梁
#             c_c = 0.85 * self.fc * self.bw * a_max
#             muc = c_c * (self.deff - a_max / 2) * self.phi
#             mus = mu - muc
#             fs = min(self.Es * self.epsilon_c_max * (c_max - self.cover) / c_max, self.fy)
#             as_c_req = mus / (fs - 0.85 * self.fc) / (self.deff - self.cover) / self.phi
#             as_t1_req = muc / self.fy / (self.deff - a_max / 2) / self.phi
#             as_t2_req = mus / self.fy / (self.deff - self.cover) / self.phi
#             as_t_req = as_t1_req + as_t2_req

#         # as_min = min(self.ps_min * self.bw * self.dp, 4 / 3 * as_t_req)
#         # as_max = self.ps_max * self.bw * self.dp

#         return as_t_req /self.bw / self.deff

#     def __ps_min(self):
#         p_min1 = 0.8 * math.sqrt(self.fc) / self.fy
#         p_min2 = 14 / self.fy
#         return max(p_min1, p_min2)

#     def __ps_max(self):
#         # 拉力鋼筋
#         if self.is_seismic_zone is True:
#             p_max1 = (self.fc + 100) / (4 * self.fy)
#             p_max2 = 0.025
#             return min(p_max1, p_max2)

#         else:
#             c = self.epsilon_c_max / (0.004 + self.epsilon_c_max)
#             a = c * self.beta1
#             return 0.85 * self.fc * a / self.fy

#     def rebar_tension_strain(self, neutral_axis_depth):
#         return self.epsilon_c_max * (self.deff - neutral_axis_depth) / neutral_axis_depth

# ----------------------------------------------------------


# rcd = Rcdesign(cfg)
# a = rcd.moment_to_rebar(mu=24.07 * 1000 * 100)
# print(f"ps_min = {rcd.ps_min : 8.4f}")
# print(f"ps_req = {a : 8.4f}")
# print(f"ps_max = {rcd.ps_max : 8.4f}")



cfg = Configuracoes()

nom_load = []
ult_load = []
nom_moment = []
ult_moment = []
eccentricity = []
phi_factor = []


c_value = np.arange(60,8,-0.5)
for c in c_value:

    es, fs = cfg.steel_strain_stress_at_depth(c, cfg.rebars[1]['coord_y'])
    esprime, fsprime = cfg.steel_strain_stress_at_depth(c, cfg.rebars[0]['coord_y'])
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

rebars = Rebars()
print(rebars.lowest)
print(rebars.total_area)