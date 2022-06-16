import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy import Point, Polygon, Ray


class RcMaterial:
    def __init__(self,fc, fy):
 
        self.fc = fc
        self.fy = fy
        self.ec_max = 0.003                     # concrete maximum strain (Constant)
        self.es_min = 0.005                     # 受拉鋼筋之淨拉應變
        self.Ec = 15000 * self.fc ** 0.5        # concrete modulus (kgf/cm2)
        self.Es = 2.04 * 10 ** 6                # steel modulus (kgf/cm2)
        self.ey = self.fy / self.Es             # steel strain at yield point
        self.beta1 = self.__beta1()

    def __beta1(self):
        if self.fc <= 280:
            beta1 = 0.85
        else:
            beta1 = max(0.85 - 0.05 * ((self.fc - 280) / 70), 0.65)
        return beta1

    def sigma_s(self, epsilon_s: float):
        # 拉力為正；壓力為負
        if epsilon_s > 0:
            return min(epsilon_s * self.Es, +self.fy)
        else:
            return max(epsilon_s * self.Es, -self.fy)

    def steel_strain_stress_at_depth(self, neutral_axis_depth, depth_of_interest):
        # 拉力為正；壓力為負
        strain = -self.ec_max/neutral_axis_depth * (neutral_axis_depth - depth_of_interest)
        stress = self.sigma_s(strain)
        return strain, stress

class Rebars:
    def __init__(self, rebars):
        self.rebars = rebars
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
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.area = self.width * self.height
    
    def get_coordinate(self):
        pass
# -----------------------------------------------------------
# 鋼筋資料
rebars = []
rebars.append({"area":3*5.067,"coord":(15,6.5)})
rebars.append({"area":0*5.067,"coord":(15,40)})
rebars.append({"area":3*5.067,"coord":(15,43.5)})
rbs = Rebars(rebars)

# -----------------------------------------------------------
# 斷面資料 unit:cm
sec = Rec(width=30, height=50)

# -----------------------------------------------------------
# 材料資料 unit:cm
mat = RcMaterial(fc = 280, fy = 4200)

# -----------------------------------------------------------
# 外力資料 unit:t-m
pu0 = 80.0
mu0 = 10.0


def forces_moments(neutral_axis_depth):
    
    a = mat.beta1 * neutral_axis_depth
    # area_sic : area of steel rebar in concrete compress area
    steel_area_in_compress_area = 0
    steel_axial_forces = 0
    steel_moments = 0
    for rebar in rbs.rebars :
        # multiplying the stress in each layer by the area in each layer.
        strain, stress = mat.steel_strain_stress_at_depth(neutral_axis_depth, rebar["coord"][1])
        steel_axial_forces += stress * rebar['area']
        steel_moments += stress * rebar['area'] * rebar["coord"][1]
        
        # 混凝土壓力區內之鋼筋面積        
        if rebar["coord"][1] < a:
            steel_area_in_compress_area += rebar['area']
    
    concrete_axial_forces = 0.85 * mat.fc * (a * sec.width - steel_area_in_compress_area)
    Pn = concrete_axial_forces - steel_axial_forces
    Mn = -concrete_axial_forces * a / 2 + steel_moments + Pn*sec.height/2
    Pu = phi * Pn
    Mu = phi * Mn
    ecc = Mn / Pn  
    return ecc, Pn, Pu, Mn, Mu

def strength_reduction_factor(epsilon_t, is_spiral = False):
    '''
    Strength reduction factor
    橫箍筋 0.65 ; 螺箍筋 0.70 
    '''
    phi = 0.7 if is_spiral == True else 0.65
    phi_factor = phi + 0.25 * (epsilon_t - mat.ey) / (0.005 - mat.ey)

    if epsilon_t >= 0.005:
        phi, classify = 0.9, "Tension-controlled"
    elif epsilon_t <= mat.ey:
        phi, classify = phi, "Compression-controlled"
    else:
        phi, classify = phi_factor, "Transition"    
    return phi, classify

def find_PMcurve_intersection(y_p0, x_m0, *curve_points):

    origin_point = Point(0, 0)
    check_point = Point(x_m0, y_p0)
    
    # 過原點直線
    ry = Ray(origin_point, check_point)
    poly = Polygon(*curve_points)
    
    # 找交點
    Segments = poly.sides
    for Segment in Segments:
        cross_point = ry.intersection(Segment)
        if len(cross_point) == 1:
            break    
    x_mi = cross_point[0].x.evalf()
    y_pi = cross_point[0].y.evalf()
    
    # 求安全係數

    dist_d = check_point.distance(origin_point)
    dist_c = cross_point[0].distance(origin_point)
    FS_ratio = ( dist_d / dist_c).evalf()
    
    return y_pi, x_mi, FS_ratio

nom_Axial_load = []
ult_Axial_load = []
nom_moment = []
ult_moment = []
eccentricity = []
phi_factor = []

# -----------------------------------------------------------
NUM = 50
c_value = np.linspace(0.0001, 1.5* sec.height, NUM)

# Pure Compression
# alfa: 橫箍筋 0.80 ; 螺箍筋 0.85 
alfa = 0.8
p0 = (0.85*mat.fc*(sec.area-rbs.total_area)+rbs.total_area*mat.fy)
pn = alfa*0.65*p0
print(f"pn = {pn/1000:.2f} tf")

PM_curve_points = []
for c in c_value:
    # 最大拉力鋼筋應變
    ept, fs = mat.steel_strain_stress_at_depth(c, sec.height-rbs.lowest)
    phi, classify =  strength_reduction_factor(ept)
    ecc, Pn, Pu, Mn, Mu = forces_moments(c)

    if ecc <= 1.5 * sec.height:
        # 轉換單位為 t-m
        nom_Axial_load.append(Pn/1000)
        ult_Axial_load.append(Pu/1000)
        nom_moment.append(Mn/100000)
        ult_moment.append(Mu/100000)
        eccentricity.append(ecc)
        phi_factor.append(phi)
        PM_curve_points.append(Point(Mn/100000, Pn/1000))

dict = {"ecc": eccentricity,
        "Pn": nom_Axial_load,
        "Pu": ult_Axial_load,
        "Mn": nom_moment,
        "Mu": ult_moment}

df = pd.DataFrame(dict)

# 求求交點
pi, mi, fs = find_PMcurve_intersection(pu0, mu0, *PM_curve_points)
print(f'P  = {pi} tf')
print(f'M  = {mi} tf-m')
print(f'FS = {fs}')

# -----------------------------------------------------------
plt.figure()

plt.plot(df["Mn"], df["Pn"], linestyle='--' , marker='.')
plt.plot(df["Mu"], df["Pu"], linestyle='-',)
plt.plot(mu0,pu0,'bx')

plt.xlim(0)
plt.title("P-M Interation Diagram")
plt.xlabel("M (tf-m)")
plt.ylabel("P (tf)")
plt.legend(["Nominal strength","Design strength"]) 
plt.grid(True)
plt.show()
