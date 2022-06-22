from sympy import Point, Polygon, Ray
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point, MultiPoint,LineString
from shapely import wkb
from shapely import affinity, ops

class RCMaterial:
    def __init__(self,fc, fy):
        # Concrete
        self.fc = fc
        self.Ec = 15000 * self.fc ** 0.5        # concrete modulus (kgf/cm2)
        self.eps_cu = 0.003                     # concrete maximum strain (Constant)
        self.beta1 = self.__beta1()             # Whitney block reduction coefficient
        
        # Steel
        self.fy = fy
        self.Es = 2.04 * 10 ** 6                # steel modulus (kgf/cm2)
        self.eps_yt = self.fy / self.Es         # steel strain at yield point
    

    def __beta1(self):
        
        if self.fc <= 280:
            return 0.85
        else:
            return max(0.85 - 0.05 * ((self.fc - 280) / 70), 0.65)

    def steel_stress(self, eps_s: float):
        if eps_s > 0:
            return min(eps_s * self.Es, +self.fy)
        else:
            return max(eps_s * self.Es, -self.fy)

    def steel_strain_stress_at_depth(self, neutral_axis_depth, depth_of_interest):
        # 拉力為負；壓力為正
        strain = self.eps_cu/neutral_axis_depth * (neutral_axis_depth - depth_of_interest)
        stress = self.steel_stress(strain)
        return strain, stress

class Rebars:
    def __init__(self, rebars):
        self.rebars = rebars
        self.lowest = self.__lowest()
        self.total_area = self.__total_area()
        self.coordinate = self.__coordinate()

    def __lowest(self):
        coord_ymin = min(rebar["coord"][1] for rebar in self.rebars)
        return coord_ymin 

    def __total_area(self):
        area = sum(rebar['area'] for rebar in self.rebars)
        return area

    def __coordinate(self):
        return [rebar["coord"] for rebar in self.rebars]   

    def rA(barSize): 
        '''
        鋼筋斷面積 cm2
        '''
        rA = 0
        if barSize == 10: rA =  0.7133
        if barSize == 13: rA =  1.267
        if barSize == 16: rA =  1.986
        if barSize == 19: rA =  2.865
        if barSize == 22: rA =  3.871
        if barSize == 25: rA =  5.067
        if barSize == 29: rA =  6.469
        if barSize == 32: rA =  8.143
        if barSize == 36: rA = 10.07
        if barSize == 39: rA = 12.19
        if barSize == 43: rA = 14.52
        if barSize == 57: rA = 25.79
        return rA


class RectangleSec:
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.area = self.width * self.height
        self.coordinate = self.__coordinate()

    def __coordinate(self):
        p1 = (-self.width/2,-self.height/2)
        p2 = (+self.width/2,-self.height/2)
        p3 = (+self.width/2,+self.height/2)
        p4 = (-self.width/2,+self.height/2)
        return [p1, p2, p3, p4]



def geometry_transform(geometry, x0, y0, angle_deg):
    geometry = affinity.translate(geometry, - x0, -y0)
    geometry = affinity.rotate(geometry, angle_deg, origin=(0, 0))
    return geometry

def computed_capacity_force(poly, neutral_axis_depth):

    # concrete_force
    # a = effective compressive block depth for ACI-318 concrete model
    a = mat.beta1 * neutral_axis_depth

    poly_xmin, poly_ymin, poly_xmax, poly_ymax = poly.bounds

    cut_line = LineString([(poly_xmin, poly_ymax - a),(poly_xmax, poly_ymax - a)])

    cut_polys = ops.split(poly, cut_line)
    cut_poly0_ymax = cut_polys.geoms[0].bounds[3]
    cut_poly1_ymax = cut_polys.geoms[1].bounds[3]
    compress_poly = cut_polys.geoms[0] if cut_poly0_ymax > cut_poly1_ymax else cut_polys.geoms[1]

    Pc = 0.85 * mat.fc * compress_poly.area        # kgf
    Mcx = Pc * -compress_poly.centroid.y           # kgf-cm
    Mcy = Pc * +compress_poly.centroid.x           # kgf-cm

    # steel_force
    # Yi  = strained fiber at level i
    # Ymax  = maximum compressive fiber
    # Y0  = “zero” strain fiber in cross-section
    Ymax = poly_ymax 
    Y0 = Ymax - neutral_axis_depth

    rbs_df["eps_s"] = mat.eps_cu / neutral_axis_depth * (rbs_df["Yri"] - Y0)
    rbs_df["sigma_s"] = [mat.steel_stress(i) for i in rbs_df["eps_s"]]


    dfc = []
    for _ in rbs_df["eps_s"]:
        if _ < 0:
            dfc.append(0)
        else:
        # 扣除壓力鋼筋範圍混凝土應力    
            dfc.append(-0.85 * mat.fc)

    rbs_df["dfc"] = dfc
    rbs_df["Ps"] = (rbs_df["sigma_s"] + rbs_df["dfc"]) * rbs_df["As"]
    rbs_df["Msx"] = - rbs_df["Ps"] * rbs_df["Yri"]
    rbs_df["Msy"] = + rbs_df["Ps"] * rbs_df["Xri"]

    Ps =  sum(rbs_df["Ps"])   # kgf
    Msx = sum(rbs_df["Msx"])  # kgf-cm
    Msy = sum(rbs_df["Msy"])  # kgf-cm

    Pn = Pc + Ps
    Mnx =  (Mcx + Msx) * math.cos(-angle_deg * math.pi/180) + \
           (Mcy + Msy) * math.sin(-angle_deg * math.pi/180)
    Mny = -(Mcx + Msx) * math.sin(-angle_deg * math.pi/180) + \
           (Mcy + Msy) * math.cos(-angle_deg * math.pi/180)

    return Pn, Mnx, Mny

def forces_moments(neutral_axis_depth):
    
    a = mat.beta1 * neutral_axis_depth
    # area_sic : area of steel rebar in concrete compress area
    steel_area_in_compress_area = 0
    steel_axial_forces = 0
    steel_moments = 0
    for rebar in rbs.rebars :
        # multiplying the stress in each layer by the area in each layer.
        steel_at_depth = rec.height/2 - rebar["coord"][1]
        strain, stress = mat.steel_strain_stress_at_depth(neutral_axis_depth, steel_at_depth)
        steel_axial_forces += stress * rebar['area']
        steel_moments += stress * rebar['area'] * steel_at_depth
        
        # 混凝土壓力區內之鋼筋面積        
        if steel_at_depth < a:
            steel_area_in_compress_area += rebar['area']
    
    concrete_axial_forces = 0.85 * mat.fc * (a * rec.width - steel_area_in_compress_area)
    Pn = concrete_axial_forces - steel_axial_forces
    Mn = -concrete_axial_forces * a / 2 + steel_moments + Pn*rec.height/2
    phi_Pn = phi * Pn
    phi_Mn = phi * Mn
    ecc = Mn / Pn  
    return ecc, Pn, phi_Pn, Mn, phi_Mn

def strength_reduction_factor(eps_t, phi_b, phi_c):
    '''
    phi_b  tension controlled   
    phi_c  compression controlled 
    '''
    eps_t = abs(eps_t)
    factor = phi_c + 0.25 * (eps_t - mat.eps_yt) / (0.005 - mat.eps_yt)
    if eps_t >= 0.005:
        phi, classify = phi_b,  "Tens."
    elif eps_t <= mat.eps_yt:
        phi, classify = phi_c,  "Comp."
    else:
        phi, classify = factor, "Tran."    
    return phi, classify

def pm_curve_intersection(y_p0, x_m0, *curve_points):

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
    demand = check_point.distance(origin_point)
    capacity = cross_point[0].distance(origin_point)
    fs = (demand/capacity).evalf()
    
    return y_pi, x_mi, fs

# -----------------------------------------------------------
# 坐標軸順時針為正
angle_deg = 180

# -----------------------------------------------------------
# 斷面資料 unit:cm
rec = RectangleSec(width=30, height=50)
rec_poly = Polygon(rec.coordinate)
rec_df = pd.DataFrame(list(rec_poly.exterior.coords), columns=['Xi', 'Yi'])


# --------------------------------------------------------
# 鋼筋資料
# 以斷面形心為原點中心 
rebars = []
rebars.append({"area":1*Rebars.rA(25),"coord":( 10,+18.5)})
rebars.append({"area":1*Rebars.rA(25),"coord":(  0,+18.5)})
rebars.append({"area":1*Rebars.rA(25),"coord":(-10,+18.5)})
rebars.append({"area":1*Rebars.rA(25),"coord":( 10,-18.5)})
rebars.append({"area":1*Rebars.rA(25),"coord":(  0,-18.5)})
rebars.append({"area":1*Rebars.rA(25),"coord":(-10,-18.5)})

rbs = Rebars(rebars)
rebarPoints = MultiPoint(rbs.coordinate)

_x   = [point.x for point in rebarPoints.geoms]
_y   = [point.y for point in rebarPoints.geoms]
_a   = [rb["area"] for rb in rbs.rebars]
rbs_df = pd.DataFrame(list(zip(_x,_y,_a)), columns=['Xi', 'Yi','As'])


# -----------------------------------------------------------
# 材料資料 unit:cm
mat = RCMaterial(fc = 280, fy = 4200)


# -----------------------------------------------------------
# capacity reduction factors
phi_a = 0.80  # axial compression 
phi_b = 0.90  # tension controlled   
phi_c = 0.65  # compression controlled 


# -----------------------------------------------------------
# plastic-centroid for the origin of the coordinate system
x_pc, y_pc = 0, 0


# -----------------------------------------------------------
# Concrete Coordinates and Cross-Sectional Properties
rec_poly = geometry_transform(rec_poly, x_pc, y_pc, angle_deg)
rec_poly_xmin, rec_poly_ymin, rec_poly_xmax, rec_poly_ymax = rec_poly.bounds
rec_df["Xri"] = rec_poly.exterior.coords.xy[0]
rec_df["Yri"] = rec_poly.exterior.coords.xy[1]
print(rec_df)

# Rebar Coordinates and Transformed Cross-Sectional Properties
rebarPoints = geometry_transform(rebarPoints, x_pc, y_pc, angle_deg)
rbs_df["Xri"] = [point.x for point in rebarPoints.geoms]
rbs_df["Yri"] = [point.y for point in rebarPoints.geoms]


# -----------------------------------------------------------
# 外力資料 unit:t-m
pu0 =  0.0
mu0 = 10.0

# Make PM curve   
NUM = 50

# 最外側拉力鋼筋距混凝土頂距離
rbs_lowest = rec_poly_ymax - rbs_df["Yri"].min() 
c_depth = rec_poly_ymax - rec_poly_ymin
c_value = np.linspace(1.0 * c_depth,0.0001, NUM)

pmCurve = []

for c in c_value:
    eps_t, fs = mat.steel_strain_stress_at_depth(c, rbs_lowest)
    phi, classify =  strength_reduction_factor(eps_t, phi_b, phi_c)
    Pn, Mnx, Mny= computed_capacity_force(rec_poly, c)
    pmCurve.append({"Pn": Pn/1000,
                    "Mnx":Mnx/100000,
                    "Mny":Mny/100000,
                    "phi_Pn": phi*Pn/1000, 
                    "phi_Mnx":phi*Mnx/100000,
                    "phi_Mny":phi*Mny/100000,
                    "NA_depth":c,
                    "eps_t":eps_t,
                    "phi": phi,})

# 全斷面受壓
# Pure Compression
# Pc = 0.85*mat.fc*(rec_poly.area-rbs.total_area) + rbs.total_area * mat.fy
# pmCurve.append({"Pn": Pc/1000, "Mnx": 0, "Mny": 0}) 

pmCurve_df = pd.DataFrame(pmCurve, columns=["Pn",
                                            "Mnx",
                                            "Mny",
                                            "phi_Pn",
                                            "phi_Mnx",
                                            "phi_Mny",
                                            "NA_depth",
                                            "eps_t",
                                            "phi"])


output = pmCurve_df.to_string(formatters={
    'eps_t': '{:,.4f}'.format,
    'Pn':  '{:,.3f}'.format,
    'Mnx': '{:,.3f}'.format,
    'Mny': '{:,.3f}'.format,
    'phi_Pn':  '{:,.3f}'.format,
    'phi_Mnx': '{:,.3f}'.format,
    'phi_Mny': '{:,.3f}'.format,
    'phi': '{:,.3f}'.format,
})
# print(output)




# Plot PM curve
plt.figure() 
plt.plot(pmCurve_df["Mnx"], pmCurve_df["Pn"],linestyle='-' , marker='.')
plt.plot(pmCurve_df["phi_Mnx"], pmCurve_df["phi_Pn"])
plt.plot(mu0,pu0,'rx')
plt.title("P-M Interation Diagram")
plt.xlabel("M (tf-m)")
plt.ylabel("P (tf)")
plt.legend(["Nominal strength","Design strength"]) 
plt.grid(True)
plt.show()



# ax = plt.figure().add_subplot(projection='3d')
# x = pmCurve_df["Mnx"]
# y = pmCurve_df["Mny"]
# z = pmCurve_df["Pn"]
# ax.plot(x, y, z, label='pm curve')
# ax.legend(["Nominal strength","Design strength"])
# ax.grid(True) 
# plt.show()


# nom_Axial_load = []
# phi_Axial_load = []
# nom_moment = []
# phi_moment = []
# eccentricity = []
# phi_factor = []
# epsilon_t = []
# neutral_axis_depth = []
# control_type =[]

# # -----------------------------------------------------------
# NUM = 20
# c_value = np.linspace(0.00001, 1.2* rec.height, NUM)


# # Pure Compression
# # Nominal axial compressive strength at zero eccentricity 
# p0 = 0.85*mat.fc*(rec.area-rbs.total_area)+rbs.total_area*mat.fy
# phi_p0 = 0.8*0.65*p0




# PM_curve_points = []
# for c in c_value:
#     # 最大拉力鋼筋應變
#     ept, fs = mat.steel_strain_stress_at_depth(c, rec.height/2-rbs.lowest)
#     
#     ecc, Pn, phi_Pn, Mn, phi_Mn = forces_moments(c)phi, classify =  strength_reduction_factor(ept)

#     nom_Axial_load.append(Pn/1000)
#     phi_Axial_load.append(phi_Pn/1000)
#     nom_moment.append(Mn/100000)
#     phi_moment.append(phi_Mn/100000)
#     eccentricity.append(ecc)
#     phi_factor.append(phi)
#     neutral_axis_depth.append(c)
#     epsilon_t.append(ept)
#     control_type.append(classify)
#     PM_curve_points.append(Point(phi_Mn/100000, phi_Pn/1000))

# dict = {"Pn": nom_Axial_load,
#         "Mn": nom_moment,
#         "phiPn": phi_Axial_load,
#         "phiMn": phi_moment,
#         "ecc": eccentricity,
#         "NA_depth": neutral_axis_depth,
#         "eps_t": epsilon_t,
#         "cont.": control_type,
#         "phi": phi_factor,}

# df = pd.DataFrame(dict)

# output = df.to_string(formatters={
#     'eps_t': '{:,.4f}'.format,
#     'Pn': '{:,.3f}'.format,
#     'Mn': '{:,.3f}'.format,
#     'NA_depth': '{:,.3f}'.format,
#     'ecc': '{:,.3f}'.format,
#     'phi': '{:,.3f}'.format,
# })
# print(output)

# 求求交點
# phi_pi, phi_mi, fs = pm_curve_intersection(pu0, mu0, *PM_curve_points)
# print(f'phi_Pn  = {phi_pi} tf')
# print(f'phi_Mn  = {phi_mi} tf-m')
# print(f'FS = {fs}')

# -----------------------------------------------------------
# plt.figure()

# plt.plot(df["Mn"], df["Pn"],linestyle='-' , marker='.')
# plt.plot(df["phiMn"], df["phiPn"], linestyle='--')
# plt.plot(mu0,pu0,'rx')

# plt.xlim(0)
# plt.title("P-M Interation Diagram")
# plt.xlabel("M (tf-m)")
# plt.ylabel("P (tf)")
# plt.legend(["Nominal strength","Design strength"]) 
# plt.grid(True)
# plt.show()

