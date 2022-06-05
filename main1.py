import math


class Configuracoes:
    def __init__(self, cover):
        self.is_seismic_zone = True
        self.epsilon_c_max = 0.003
        self.epsilon_s_min = 0.005  # 受拉鋼筋之淨拉應變
        self.cover = cover
        self.max_rebar_tension_strain = 0.004
        self.phi = 0.9


class Rectshape:
    def __init__(self, width, depth):
        self.width = width
        self.depth = depth
        self.area = self.width * self.depth


class Material:
    def __init__(self, fc, fy):
        self.fc = fc
        self.fy = fy
        self.Es = 2.04 * 10 ** 6  # steel modulus (kgf/cm2)
        self.Ec = 15000 * self.fc ** 0.5
        self.epsilon_y = self.fy / self.Es  # steel strain at yield point
        self.beta1 = self.__beta1()

    # Beta1 (stress block multiplier)
    def __beta1(self):
        if self.fc <= 280:
            beta1 = 0.85
        else:
            beta1 = max(0.85 - 0.05 * ((self.fc - 280) / 70), 0.65)
        return beta1

    # Strength reduction factor
    def __strength_factor(self, es):
        phi_factor = 0.65 + 0.25 * (es - self.epsilon_y) / (0.005 - self.epsilon_y)
        if phi_factor >= 0.9:
            phi, classify = 0.9, "Tension-controlled"
        elif phi_factor <= 0.65:
            phi, classify = 0.65, "Compression-controlled"
        else:
            phi, classify = phi_factor, "Transition"
        return phi, classify


class Rcdesign:
    def __init__(self, mat: Material, rec: Rectshape, cfg: Configuracoes):
        """
        矩形梁鋼筋設計
        :type mat: object
        """
        self.mat = mat
        self.rec = rec
        self.cfg = cfg
        self.fc = mat.fc
        self.fy = mat.fy
        self.Es = mat.Es
        self.beta1 = mat.beta1
        self.bw = rec.width
        self.depth = rec.depth
        self.phi = cfg.phi
        self.cover_compres = cfg.cover
        self.cover_tension = cfg.cover
        self.epsilon_c_max = cfg.epsilon_c_max
        self.epsilon_s_min = cfg.epsilon_s_min
        self.is_seismic_zone = cfg.is_seismic_zone
        self.dp = self.depth - self.cover_tension
        self.ps_min = self.__ps_min()
        self.ps_max = self.__ps_max()
        self.as_min = self.ps_min * self.bw * self.dp
        self.as_max = self.ps_max * self.bw * self.dp

    def moment_to_rebar(self, mu: float):
        """
        :param mu: 係數化彎矩 kgf-cm
        :return:
        """
        # a : the depth of the compression block
        a = self.dp - (self.dp ** 2 - 2 * abs(mu) / (0.85 * self.fc * self.phi * self.bw)) ** 0.5
        c_max = self.epsilon_c_max / (self.epsilon_c_max + self.epsilon_s_min) * self.dp
        a_max = self.beta1 * c_max

        if a <= a_max:
            # 單筋梁
            as_t_req = mu / (self.phi * self.fy * (self.dp - a / 2))
            as_c_req = 0

        else:
            # 雙筋梁
            c_c = 0.85 * self.fc * self.bw * a_max
            muc = c_c * (self.dp - a_max / 2) * self.phi
            mus = mu - muc
            fs = min(self.Es * self.epsilon_c_max * (c_max - self.cover_compres) / c_max, self.fy)
            as_c_req = mus / (fs - 0.85 * self.fc) / (self.dp - self.cover_compres) / self.phi
            as_t1_req = muc / self.fy / (self.dp - a_max / 2) / self.phi
            as_t2_req = mus / self.fy / (self.dp - self.cover_compres) / self.phi
            as_t_req = as_t1_req + as_t2_req

        # as_min = min(self.ps_min * self.bw * self.dp, 4 / 3 * as_t_req)
        # as_max = self.ps_max * self.bw * self.dp

        return as_t_req / self.bw / self.dp

    def __ps_min(self):
        p_min1 = 0.8 * math.sqrt(self.fc) / self.fy
        p_min2 = 14 / self.fy
        return max(p_min1, p_min2)

    def __ps_max(self):
        # 拉力鋼筋
        if self.is_seismic_zone is True:
            p_max1 = (self.fc + 100) / (4 * self.fy)
            p_max2 = 0.025
            return min(p_max1, p_max2)

        else:
            c = self.epsilon_c_max / (0.004 + self.epsilon_c_max)
            a = c * self.beta1
            return 0.85 * self.fc * a / self.fy

    def rebar_tension_strain(self, neutral_axis_depth):
        return self.epsilon_c_max * (self.dp - neutral_axis_depth) / neutral_axis_depth


# ----------------------------------------------------------

mat = Material(fc=280, fy=4200)
rec = Rectshape(width=35, depth=70)
cfg = Configuracoes(cover=7.2)
rcd = Rcdesign(mat, rec, cfg)
a = rcd.moment_to_rebar(mu=24.07 * 1000 * 100)
print(f"ps_min = {rcd.ps_min : 8.4f}")
print(f"ps_req = {a : 8.4f}")
print(f"ps_max = {rcd.ps_max : 8.4f}")
