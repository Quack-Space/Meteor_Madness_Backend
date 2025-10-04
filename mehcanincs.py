#example_var_name

#from poliastro import core as poliore
import poliastro
import numpy as np
import utm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

################
### CONSTANTS ###

grav_constant = 6.674 * 10**-11         #N·m²/kg²

R_Earth  = 6378137.0;                   #[m] medium equatorial radius
earth_mass = 5.972 * 10**24             #[kg]

mu_Earth = grav_constant * earth_mass   #[m^3/s^2] Earth gravitation constant #398600441800000.00

# Sun constants for Sun-Heliocentric Reference (SHRC)
sun_mass = 1.98847e30                    # [kg]
mu_Sun = grav_constant * sun_mass        # [m^3/s^2] Sun gravitational parameter ~1.327e20

J2 = 1.08262668e-3;                     # [-] oblateness coefficient
J3 = -2.53215306e-6;                    # [-] oblateness coefficient
J4 = -1.61098761e-6;                    # [-] oblateness coefficient

################

#Convert latitude and longitude to UTM coordinates.
def latlon_2_utm(lat, lon):
    utm_coords = utm.from_latlon(np.array(lat), np.array(lon))
    return utm_coords

#Convert UTM coordinates to latitude and longitude.
def utm_2_latlon(easting, northing, zone_number, zone_letter):
    lat, lon = utm.to_latlon(easting, northing, zone_number, zone_letter)
    return lat, lon

#Convert latitude and longitude to ECEF coordinates.
def latlon_2_ecef(lat, lon, alt):
    lat = np.radians(lat)
    lon = np.radians(lon)
    a = 6378137.0
    e = 8.1819190842622e-2
    N = a / np.sqrt(1 - e**2 * np.sin(lat)**2)
    x = (N + alt) * np.cos(lat) * np.cos(lon)
    y = (N + alt) * np.cos(lat) * np.sin(lon)
    z = (N * (1 - e**2) + alt) * np.sin(lat)
    return x, y, z

#Convert ECEF coordinates to latitude and longitude.
def ecef_2_latlon(x, y, z):
    a = 6378137.0
    e = 8.1819190842622e-2
    b = np.sqrt(a**2 * (1 - e**2))
    ep = np.sqrt((a**2 - b**2) / b**2)
    p = np.sqrt(x**2 + y**2)
    th = np.arctan2(a * z, b * p)
    lon = np.arctan2(y, x)
    lat = np.arctan2((z + ep**2 * b * np.sin(th)**3), (p - e**2 * a * np.cos(th)**3))
    N = a / np.sqrt(1 - e**2 * np.sin(lat)**2)
    alt = p / np.cos(lat) - N
    lat = np.degrees(lat)
    lon = np.degrees(lon)
    return lat, lon, alt

#Convert ECEF coordinates to ECI coordinates.
def ecef_2_eci(x, y, z, gst):
    gst = np.radians(gst)
    x_eci = x * np.cos(gst) - y * np.sin(gst)
    y_eci = x * np.sin(gst) + y * np.cos(gst)
    z_eci = z
    return x_eci, y_eci, z_eci

#Convert ECI coordinates to ECEF coordinates.
def eci_2_ecef(x, y, z, gst):
    gst = np.radians(gst)
    x_ecef = x * np.cos(gst) + y * np.sin(gst)
    y_ecef = -x * np.sin(gst) + y * np.cos(gst)
    z_ecef = z
    return x_ecef, y_ecef, z_ecef
################

# Keplerian -> Cartesian (SHRF)
#Convert Keplerian orbital elements to SHRF position vector (meters).

#Inputs:
#- a: semi-major axis [m]
#- e: eccentricity (0<=e<1)
#- i_deg: inclination [deg]
#- raan_deg: right ascension of ascending node Omega [deg]
#- argp_deg: argument of perigee omega [deg]
#- nu_deg: true anomaly [deg]

#Returns: (x, y, z) in SHRF frame [m]
def keplerian_to_shrf(a, e, i_deg, raan_deg, argp_deg, nu_deg):
    if e >= 1.0:
        raise ValueError("This function currently supports only elliptical orbits (e < 1).")

    #deg 2 rad
    i = np.radians(i_deg)
    raan = np.radians(raan_deg)
    argp = np.radians(argp_deg)
    nu = np.radians(nu_deg)

    #calculating distance in orbital plane
    r = a * (1 - e**2) / (1 + e * np.cos(nu))

    #position in perifocal coordinates
    x_pf = r * np.cos(nu)
    y_pf = r * np.sin(nu)
    z_pf = 0.0

    #Rotation matrices
    #R = Rz(raan) * Rx(i) * Rz(argp)
    ca = float(np.cos(raan)); sa = float(np.sin(raan))
    ci = float(np.cos(i));    si = float(np.sin(i))
    cw = float(np.cos(argp)); sw = float(np.sin(argp))

    # Combined rotation matrix from perifocal to ECI
    R = np.array([
        [ca * cw - sa * sw * ci, -ca * sw - sa * cw * ci, sa * si],
        [sa * cw + ca * sw * ci, -sa * sw + ca * cw * ci, -ca * si],
        [sw * si               , cw * si                , ci]
        ])

    x_eci = R[0, 0] * x_pf + R[0, 1] * y_pf + R[0, 2] * z_pf
    y_eci = R[1, 0] * x_pf + R[1, 1] * y_pf + R[1, 2] * z_pf
    z_eci = R[2, 0] * x_pf + R[2, 1] * y_pf + R[2, 2] * z_pf
    return float(x_eci), float(y_eci), float(z_eci)


#Generate a list of (x,y,z) tuples representing points on an elliptical orbit
#in the Sun-Heliocentric Reference Frame (SHRF) defined by Keplerian elements.
#You can specify either a delta in true anomaly (degrees) with `delta_nu_deg`,
#or an approximate arc-length spacing in meters with `spacing`. Alternatively,
#request `num_points` equally spaced in true anomaly.

#Priority of options: if `delta_nu_deg` is provided it is used; else if
#`num_points` provided produce that many points; else default to 361 points
#(1 deg steps).

#Returns: list of (x,y,z) tuples in SHRF frame (meters)

def generate_ellipse_points_shrf(a, e, i_deg=0.0, raan_deg=0.0, argp_deg=0.0,start_nu_deg=0.0, delta_nu_deg=None,num_points=None):
    if e >= 1.0:
        raise ValueError("Only elliptical orbits (e < 1) are supported.")
    start = start_nu_deg % 360.0

    if delta_nu_deg is not None:
        if delta_nu_deg <= 0:
            raise ValueError("delta_nu_deg must be positive")
        nus = np.arange(start, start + 360.0, delta_nu_deg)
        points = [keplerian_to_shrf(a, e, i_deg, raan_deg, argp_deg, float(nu))
                  for nu in nus]
        return points

    # fallback: num_points or default
    if num_points is not None:
        if num_points <= 0:
            raise ValueError("num_points must be positive")
        nus = np.linspace(start, start + 360.0, num_points, endpoint=False)
    else:
        nus = np.linspace(start, start + 360.0, 361, endpoint=False)

    points = [keplerian_to_shrf(a, e, i_deg, raan_deg, argp_deg, float(nu))
              for nu in nus]
    return points



def plot_points_3d(points, title="Orbit points", outpath="orbit_plot.png"):
        pts_arr = np.array(points)
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(pts_arr[:,0], pts_arr[:,1], pts_arr[:,2], lw=1)
        ax.scatter(pts_arr[0,0], pts_arr[0,1], pts_arr[0,2], color='red', s=20, label='start')
        ax.set_xlabel('X [m]')
        ax.set_ylabel('Y [m]')
        ax.set_zlabel('Z [m]')
        ax.set_title(title)
        ax.legend()
        fig.savefig(outpath, dpi=200)
        plt.close(fig)


if __name__ == "__main__":
    # quick example / smoke test
    # Example: low Earth orbit ~7000 km semi-major, small eccentricity
    a = 149.6e6
    e = 0.01671123
    i = 1.578690
    raan = 174.9
    argp = 288.1

    pts = generate_ellipse_points_shrf(a, e, i_deg=i, raan_deg=raan, argp_deg=argp, num_points=1000) #delta_nu_deg=30.0)
    print(f"Generated {len(pts)} points.")

    out_file = "orbit_shrf.png"
    plot_points_3d(pts, title="SHRF Orbit", outpath=out_file)
    print(f"Saved plot to {out_file}")








##############################################################################################################################







from poliastro import core as policore
import numpy as np
import math


#J2 perturbations
#state comprehends the six component state vector [x,y,z,vx,vy,vz]
def J2(t0, state, k, J2, r):
    J2pert = policore.perturbations.J2_perturbation(t0, state, k, J2, r)
    return J2pert


#atm drag acceleration
def atm_drag(t0, state, k, R, C_D, A_over_m, rho):
    adrag_acc = policore.perturbations.atmospheric_drag(t0, state, k, R, C_D, A_over_m, rho)
    return adrag_acc

#radiation pressure acceleration
def rad_pressure(t0, state, k, r, C_R, A_over_m, Wdivc_s, star):
    rad_p = policore.perturbations.radiation_pressure(t0, state, k, R, C_R, A_over_m, Wdivc_s, star)
    return rad_p

#final velocity (V entry atm)
def final_v(v0, angle, h, m, area, cd, rho, g=9.81):
    k = 0.5 * rho * cd * area
    theta = math.radians(angle)

    v0x = v0 * math.cos(theta)
    v0y = -v0 * math.sin(theta)  # verso il basso (negativo)

    exp_kh = math.exp(-k * h / m)
    exp_2kh = math.exp(-2 * k * h / m)

    vfx = v0x * exp_kh
    vfy_squared = (m * g / k) * (1 - exp_2kh) + v0y**2 * exp_2kh
    vfy = math.sqrt(vfy_squared)
    vf = math.sqrt(vfx**2 + vfy**2)
    return vf

#Ek from mass/rdot
def kinetic_energy(m, v, vf):
    ek = 0.5*m*(v**2)
    ekf = 0.5*m*(vf**2)
    return ek, ekf


#crater dimensions  
# Constants for target materials
target_par = {
    "rock": {"K1": 0.22, "mu": 0.55, "K2": 0.2, "nu": 0.4, "density": 2700, "strength": 1e7},  # strength in Pa (~10 MPa)
    "soil": {"K1": 0.24, "mu": 0.41, "K2": 0.55, "nu": 0.33, "density": 1800, "strength": 1e6},
    "ice": {"K1": 0.19, "mu": 0.55, "K2": 0.24, "nu": 0.4, "density": 900, "strength": 5e6},
    }

#here m, v and d (diameter) are the initial values (before the atm transit)
def atm_transit(m, v, d, density_impactor,
                        impact_angle,
                        sigma=1e-4,  # ablation coeff (1/mass per area)
                        kappa=1e-5,  # velocity loss coeff
                        atm_column_density=1.0e4):  # kg/m^2 at sea level vertical
    
    #Simple atmospheric ablation and velocity loss
    #mass and velocity exponentially decay with atmospheric column density traversed.
    
    theta = math.radians(impact_angle)
    if math.sin(theta) < 1e-3:
        sin_theta = 1e-3  # avoid division by zero for very shallow angles
    else:
        sin_theta = math.sin(theta)

    X_theta = atm_column_density / sin_theta  # effective atmospheric column along trajectory

    mf = m * math.exp(-sigma * X_theta)
    vf = v * math.exp(-kappa * X_theta)

    # Diameter scales with cube root of mass assuming constant density
    af = d * (mf / m) ** (1/3)

    return mf, vf, af

def crater_with_angle_correction(Dt, impact_angle):
    theta = math.radians(impact_angle)
    # Angle correction exponent can be tuned
    angle_factor = (math.sin(theta))**(1/3)
    return Dt * angle_factor

def compute_transient_crater_diameter(mf, vf, af, density_impactor,
                                      target_density, g, target_strength,
                                      K1, mu, K2, nu):
   
    #Compute transient crater diameter Dt using π-scaling (strength + gravity).
    
    pi2 = g * af / (vf**2)
    pi3 = target_strength / (target_density * vf**2)

    piD = K1 * (pi2 ** (-mu)) + K2 * (pi3 ** (-nu))

    Dt = (vf**2 / g) * piD
    return Dt

def collapse_correction(Dt, simple_to_complex_D, kd_simple=1.2, kd_complex_exponent=0.3):
    if Dt <= simple_to_complex_D:
        return Dt * kd_simple
    else:
        return Dt * (kd_simple * (Dt / simple_to_complex_D) ** kd_complex_exponent)


def depth_from_diameter(D, depth_ratio_simple=0.2, depth_ratio_complex=0.1):
    if D <= 3000:
        return depth_ratio_simple * D
    else:
        return depth_ratio_complex * D
    
def crater_dimensions_advanced(m, v, d,
                               density_impactor, impact_angle,
                               target_type="rock", g=9.81,
                               simple_to_complex_D=3000,
                               sigma=1e-4, kappa=1e-5, atm_column_density=1.0e4,
                               depth_ratio_simple=0.2, depth_ratio_complex=0.1):
    
    #Full function with material choice and ablation included.
    

    if target_type not in target_par:
        raise ValueError(f"Unknown target type '{target_type}', valid options: {list(target_par.keys())}")

    params = target_par[target_type]

    # Atmospheric transit with ablation and velocity loss
    m_f, v_f, a_f = atm_transit(m, v, d,
                                        density_impactor, impact_angle,
                                        sigma=sigma, kappa=kappa, atm_column_density=atm_column_density)
    # Transient crater diameter
    D_t = compute_transient_crater_diameter(m_f, v_f, a_f, density_impactor,
                                            params["density"], g, params["strength"],
                                            params["K1"], params["mu"], params["K2"], params["nu"])

    # Angle correction
    Dt_corr = crater_with_angle_correction(D_t, impact_angle)

    # Collapse correction (simple/complex)
    D_final = collapse_correction(Dt_corr, simple_to_complex_D)

    # Depth estimate
    depth = depth_from_diameter(D_final, depth_ratio_simple, depth_ratio_complex)

    return D_final, depth
