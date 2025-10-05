import poliastro.core as policore
import numpy as np
import math
import utm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import kingery_bulmash as kb

################
### CONSTANTS ###

grav_constant = 6.674e-11         #N·m²/kg²

R_Earth  = 6378137.0;                   #[m] medium equatorial radius

earth_mass = 5.972 * 10**24             #[kg]
mu_Earth = grav_constant * earth_mass   #[m^3/s^2] Earth gravitation constant #398600441800000.00

# Sun constants for Sun-Heliocentric Reference Frame(SHRF)
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

#Static orbuit determination SOD
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
    if e >= 1.0 and e < 0.0:
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

    #Combined rotation matrix from perifocal to SHRF
    R = np.array([
        [ca * cw - sa * sw * ci, -ca * sw - sa * cw * ci, sa * si],
        [sa * cw + ca * sw * ci, -sa * sw + ca * cw * ci, -ca * si],
        [sw * si               , cw * si                , ci]
        ])

    x_shrf = R[0, 0] * x_pf + R[0, 1] * y_pf + R[0, 2] * z_pf
    y_shrf = R[1, 0] * x_pf + R[1, 1] * y_pf + R[1, 2] * z_pf
    z_shrf = R[2, 0] * x_pf + R[2, 1] * y_pf + R[2, 2] * z_pf
    return round(float(x_shrf),1), round(float(y_shrf),1), round(float(z_shrf),1)


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
    if e >= 1.0 and e < 0.0:
        raise ValueError("Only elliptical orbits (e < 1) are supported.")
    start = start_nu_deg % 360.0

    if delta_nu_deg is not None:
        if delta_nu_deg <= 0:
            raise ValueError("delta_nu_deg must be positive")
        nus = np.arange(start, start + 360.0, delta_nu_deg)
        points = [keplerian_to_shrf(a, e, i_deg, raan_deg, argp_deg, float(nu))
                  for nu in nus]
        return points

    if num_points is not None:
        if num_points <= 0:
            raise ValueError("num_points must be positive")
        nus = np.linspace(start, start + 360.0, num_points, endpoint=False)
    else:
        nus = np.linspace(start, start + 360.0, 361, endpoint=False)

    points = [keplerian_to_shrf(a, e, i_deg, raan_deg, argp_deg, float(nu))
              for nu in nus]
    return points






#DYNAMIC ORBIT DETERMINATION DOD
#Keplerian -> Cartesian (SHRF)
#Convert Keplerian orbital elements to SHRF position vector (meters).

#Inputs:
#- a: semi-major axis [m]
#- e: eccentricity (0<=e<1)
#- i_deg: inclination [deg]
#- raan_deg: right ascension of ascending node Omega [deg]
#- argp_deg: argument of perigee omega [deg]
#- nu_deg: true anomaly [deg]
#Compute position and velocity in SHRF (heliocentric inertial) frame for
#a given set of Keplerian elements and true anomaly.

#Returns: (x, y, z),(vx, vy, vz) in SHRF frame [m] [m/s]
def keplerian_to_timed_shrf(a, e, i_deg, raan_deg, argp_deg, nu_deg, mu = mu_Sun):
    if e >= 1.0 and e < 0.0:
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

    #semi-latus rectum
    p = a * (1 - e**2)

    #specific angular momentum
    h = np.sqrt(mu * p)

    #velocity in perifocal
    vx_pf = -mu / h * np.sin(nu)
    vy_pf = mu / h * (e + np.cos(nu))
    vz_pf = 0.0

    #Rotation matrices
    #R = Rz(raan) * Rx(i) * Rz(argp)
    ca = float(np.cos(raan)); sa = float(np.sin(raan))
    ci = float(np.cos(i));    si = float(np.sin(i))
    cw = float(np.cos(argp)); sw = float(np.sin(argp))

    #Combined rotation matrix from perifocal to SHRF
    R = np.array([
        [ca * cw - sa * sw * ci, -ca * sw - sa * cw * ci, sa * si],
        [sa * cw + ca * sw * ci, -sa * sw + ca * cw * ci, -ca * si],
        [sw * si               , cw * si                , ci]
        ])

    x = R[0, 0] * x_pf + R[0, 1] * y_pf + R[0, 2] * z_pf
    y = R[1, 0] * x_pf + R[1, 1] * y_pf + R[1, 2] * z_pf
    z = R[2, 0] * x_pf + R[2, 1] * y_pf + R[2, 2] * z_pf

    vx = R[0, 0] * vx_pf + R[0, 1] * vy_pf + R[0, 2] * vz_pf
    vy = R[1, 0] * vx_pf + R[1, 1] * vy_pf + R[1, 2] * vz_pf
    vz = R[2, 0] * vx_pf + R[2, 1] * vy_pf + R[2, 2] * vz_pf

    return (round(float(x), 1), round(float(y), 1), round(float(z), 1)), (round(float(vx), 1), round(float(vy), 1), round(float(vz), 1))


#Compute time since periapsis passage (seconds) for a given true anomaly nu
#using Kepler's equation. Returns time in seconds (>=0 and < orbital period).
def true_anomaly_to_time_since_periapsis(a, e, nu_deg, mu=mu_Sun):
    if e >= 1.0 and e < 0.0:
        raise ValueError("Only elliptical orbits (e < 1) are supported.")

    #deg 2 rad
    nu = np.radians(nu_deg)
    #eccentric anomaly E from true anomaly
    if abs(e) < 1e-12:
        E = nu
    else:
        tan_halfE = np.sqrt((1 - e) / (1 + e)) * np.tan(nu / 2.0)
        E = 2.0 * np.arctan(tan_halfE)
        if E < 0:
            E = E + 2.0 * np.pi

    #mean anomaly
    M = E - e * np.sin(E)

    #mean motion
    n = np.sqrt(mu / (a ** 3))

    #time since periapsis
    t = M / n
    #orbital period
    T = 2.0 * np.pi / n
    #normalize to [0, T)
    t = t % T
    return float(t)


#Generate positions, velocities, and timestamps (seconds since epoch or periapsis)
#for points along an elliptical heliocentric orbit.

#Returns a tuple (positions, velocities, times):
#- positions: list of (x,y,z) in meters (SHRF inertial)
#- velocities: list of (vx,vy,vz) in m/s
#- times: list of timestamps in seconds. If `epoch` is None, times are
#seconds since periapsis (>=0). If `epoch` is a float, it is treated
#as seconds since some reference and `times` = epoch + t_since_peri.
#If `epoch` is a datetime.datetime, it is converted to POSIX timestamp.

#Generate a list of (x,y,z) tuples representing points on an elliptical orbit
#in the Sun-Heliocentric Reference Frame (SHRF) defined by Keplerian elements.
#You can specify either a delta in true anomaly (degrees) with `delta_nu_deg`,
#or an approximate arc-length spacing in meters with `spacing`. Alternatively,
#request `num_points` equally spaced in true anomaly.

#Priority of options: if `delta_nu_deg` is provided it is used; else if
#`num_points` provided produce that many points; else default to 361 points
#(1 deg steps).

#Returns: list of (x,y,z) tuples in SHRF frame (meters)

def generate_ellipse_timed_points_shrf(a, e, i_deg=0.0, raan_deg=0.0, argp_deg=0.0,start_nu_deg=0.0, mu = mu_Sun, delta_nu_deg=None, num_points=None, epoch=None):
    if e >= 1.0 and e < 0.0:
        raise ValueError("Only elliptical orbits (e < 1) are supported.")
    start = start_nu_deg % 360.0
    positions = []
    velocities = []
    times_rel = []

    if delta_nu_deg is not None:
        if delta_nu_deg <= 0:
            raise ValueError("delta_nu_deg must be positive")
        nus = np.arange(start, start + 360.0, delta_nu_deg)
        Datas = [(keplerian_to_timed_shrf(a, e, i_deg, raan_deg, argp_deg, float(nu)), true_anomaly_to_time_since_periapsis(a, e, float(nu), mu)) for nu in nus]

        [positions.append(_[0][0]) for _ in Datas]
        [velocities.append(_[0][1]) for _ in Datas]
        [times_rel.append(round(_[1],1)) for _ in Datas]
        
        #convert epoch if provided
        if epoch is None:
            times = times_rel
        else:
            times = [round((float(epoch) + t), 1) for t in times_rel]

        return positions, velocities, times
    elif num_points is not None:
        if num_points <= 0:
            raise ValueError("num_points must be positive")
        nus = np.linspace(start, start + 360.0, num_points, endpoint=False)
    else:
        nus = np.linspace(start, start + 360.0, 361, endpoint=False)

    Datas = [(keplerian_to_timed_shrf(a, e, i_deg, raan_deg, argp_deg, float(nu)), true_anomaly_to_time_since_periapsis(a, e, float(nu), mu)) for nu in nus]

    [positions.append(_[0][0]) for _ in Datas]
    [velocities.append(_[0][1]) for _ in Datas]
    [times_rel.append(round(_[1],1)) for _ in Datas]
    
    #convert epoch if provided
    if epoch is None:
        times = times_rel
    else:
        times = [round((float(epoch) + t), 1) for t in times_rel]

    return positions, velocities, times
    


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



###EXAMPLE USAGE OF STATIC AND DYNAMIC ORBIT POINTS CALCULATION
"""
if __name__ == "__main__":
    # quick example / smoke test
    a = 149.6e9
    e = 0.017
    i = 1.58
    raan = 174.9
    argp = 288.1

    a2 = 227.94e9
    e2 = 0.093
    i2 = 1.85
    raan2 = 149.6
    argp2 = 286.5

    pos, vel, time = generate_ellipse_timed_points_shrf(a, e, i_deg=i, raan_deg=raan, argp_deg=argp, epoch=1759613269, num_points=1000)#, delta_nu_deg=30.0, num_points=10
    print(f"Generated {len(pos)} points.")
    
    print("pos")
    for p in pos:
        print (p)
    print("vel")
    for v in vel:
        print (v)
    print("time")
    for t in time:
        print (t)
    out_file = "orbit_timed_shrf.png"
    plot_points_3d(pos, title="SHRF Orbit", outpath=out_file)
    print(f"Saved plot to {out_file}")

#Static plotting 1
    pts = generate_ellipse_points_shrf(a, e, i_deg=i, raan_deg=raan, argp_deg=argp, num_points=1000) #delta_nu_deg=30.0)
    print(f"Generated {len(pts)} points.")
#    print("pts1")
#    for p in pts:
#        print (p)
    pts_arr = np.array(pts)
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(pts_arr[:,0], pts_arr[:,1], pts_arr[:,2], lw=1)
    ax.scatter(pts_arr[0,0], pts_arr[0,1], pts_arr[0,2], color='red', s=20, label='start')

#Static plotting 2
    pts = generate_ellipse_points_shrf(a2, e2, i_deg=i2, raan_deg=raan2, argp_deg=argp2, num_points=1000) #delta_nu_deg=30.0)
    print(f"Generated {len(pts)} points.")

#    print("pts2")
#    for p in pts:
#        print (p)

    pts_arr = np.array(pts)
    ax.plot(pts_arr[:,0], pts_arr[:,1], pts_arr[:,2], lw=1)
    ax.scatter(pts_arr[0,0], pts_arr[0,1], pts_arr[0,2], color='red', s=20, label='start')

#Def before saving
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')
    ax.set_title("SHRF static Orbit")
    ax.legend()

    fig.savefig("orbit_static_shrf.png", dpi=200)

    print(f"Saved plot to orbit_static_shrf.png")
"""


##############################################################################################################################


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
    rad_p = policore.perturbations.radiation_pressure(t0, state, k, r, C_R, A_over_m, Wdivc_s, star)
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
"""
target_par = {
    "rock": {"K1": 0.22, "mu": 0.55, "K2": 0.2, "nu": 0.4, "density": 2700, "strength": 1e7},  # strength in Pa (~10 MPa)
    "soil": {"K1": 0.24, "mu": 0.41, "K2": 0.55, "nu": 0.33, "density": 1800, "strength": 1e6},
    "ice": {"K1": 0.19, "mu": 0.55, "K2": 0.24, "nu": 0.4, "density": 900, "strength": 5e6},
    }
"""
target_par = {
    "Unconsolidated sediments": {"K1": 0.24, "mu": 0.41, "K2": 0.55, "nu": 0.33, "density": 1800, "strength": 1e6},
    "Siliciclastic sedimentary rocks": {"K1": 0.22, "mu": 0.55, "K2": 0.2, "nu": 0.4, "density": 2500, "strength": 1e7},
    "Mixed sedimentary rocks": {"K1": 0.22, "mu": 0.55, "K2": 0.2, "nu": 0.4, "density": 2600, "strength": 1e7},
    "Carbonate sedimentary rocks": {"K1": 0.23, "mu": 0.55, "K2": 0.21, "nu": 0.4, "density": 2700, "strength": 1e7},
    "Metamorphics": {"K1": 0.25, "mu": 0.6, "K2": 0.22, "nu": 0.35, "density": 2800, "strength": 2e7},
    "Acid/Intermediate/Basic plutonic rocks": {"K1": 0.26, "mu": 0.6, "K2": 0.23, "nu": 0.35, "density": 2800, "strength": 2e7},
    "Acid/Intermediate/Basic volcanic rocks": {"K1": 0.25, "mu": 0.55, "K2": 0.22, "nu": 0.35, "density": 2700, "strength": 1.5e7},
    "Pyroclastics": {"K1": 0.24, "mu": 0.5, "K2": 0.21, "nu": 0.35, "density": 2300, "strength": 1e7},
    "Evaporites": {"K1": 0.2, "mu": 0.4, "K2": 0.18, "nu": 0.33, "density": 2200, "strength": 5e6},
    "Ice and Glaciers": {"K1": 0.19, "mu": 0.55, "K2": 0.24, "nu": 0.4, "density": 900, "strength": 5e6},
    "Water Bodies": {"K1": 0.0, "mu": 0.0, "K2": 0.0, "nu": 0.5, "density": 1000, "strength": 0}
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
                               target_type="Mixed sedimentary rocks", g=9.81,
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



##############################################################################################################################



#DAMAGE CALCULATOR
# Requires: pip install kingery-bulmash

# constants
PSI_TO_KPA = 6.89476
J_PER_KG_TNT = 4.184e6  # J per kg TNT

#Takes kinetic energy in Joules and terrain seismic efficiency (eta)
#ONLY IF GROUND IMPACT, NOT WATER
def get_seismic_equivalent(kinetic_energy_joules, eta):
        kinetic_energy_transfered = eta * kinetic_energy_joules
        magnitude = ((2/3) * (math.log10(kinetic_energy_transfered) - 4.8))/1.5
        return magnitude


def energy_to_kgTNT(E_joules: float) -> float:
    """Convert blast energy (J) to kg TNT."""
    return E_joules / J_PER_KG_TNT


def incident_pressure_kpa(neq_kg: float, distance_m: float) -> float:
    """
    Return incident overpressure (kPa) from kingery_bulmash for a hemispherical surface burst.
    - neq_kg: TNT equivalent mass in kg
    - distance_m: distance from burst in meters
    """
    # safe=True will raise ValueError if out-of-range; we rely on exceptions to handle invalid distances
    res = kb.Blast_Parameters(unit_system=kb.Units.METRIC, neq=neq_kg, distance=distance_m, safe=False)
    return float(res.incident_pressure)  # kPa


def find_radius_for_overpressure(neq_kg: float, target_psi: float, r_lo=0.1, r_hi=None, tol=1e-2, max_iter=60):
    """
    Find radius R (m) such that incident_pressure(R) ~= target_psi, using bisection.
    - neq_kg: TNT mass in kg
    - target_psi: overpressure target in psi
    - r_lo: initial lower bound (m)
    - r_hi: initial upper bound (m); if None it's estimated from scaled-distance bounds
    - tol: convergence tolerance (meters)
    """

    target_kpa = target_psi * PSI_TO_KPA

    # estimate sensible upper bound if not provided using Kingery-Bulmash valid scaled-distance (Z_max ~40 m/kg^(1/3))
    if r_hi is None:
        # W^(1/3) * Z_max => rough upper bound for distance where KB is still valid
        Z_MAX = 40.0  # m / kg^(1/3) (KB valid up to ~40)
        r_hi = Z_MAX * (neq_kg ** (1.0/3.0))

    # Make sure bounds produce function values that bracket the root (pressure decreases with distance)
    try:
        p_lo = incident_pressure_kpa(neq_kg, r_lo)
    except Exception:
        p_lo = float('inf')  # treat as very large if invalid (very near-field)
    p_hi = incident_pressure_kpa(neq_kg, r_hi)

    # If even at r_lo pressure is below the target, no root in (r_lo,r_hi) — try even smaller r_lo
    if p_lo < target_kpa:
        # try shrinking r_lo a bit
        r_lo = r_lo * 0.1
        try:
            p_lo = incident_pressure_kpa(neq_kg, r_lo)
        except Exception:
            p_lo = float('inf')

    # If at the estimated r_hi pressure is still > target, extend r_hi (rare for huge W or tiny target)
    iter_extend = 0
    while p_hi > target_kpa and iter_extend < 10:
        r_hi *= 2.0
        p_hi = incident_pressure_kpa(neq_kg, r_hi)
        iter_extend += 1

    if not (p_lo >= target_kpa and p_hi <= target_kpa):
        # If we can't bracket the target, return None to indicate failure
        return None

    # Bisection loop
    lo, hi = r_lo, r_hi
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        try:
            p_mid = incident_pressure_kpa(neq_kg, mid)
        except Exception:
            # if res returns None/invalid in extreme near-field, move lo inward
            lo = max(lo * 0.1, 1e-6)
            continue

        if abs(p_mid - target_kpa) <= 1e-2:  # ~0.01 kPa tolerance
            return mid
        # pressure monotonically decreases with distance: if p_mid > target -> move lo up
        if p_mid > target_kpa:
            lo = mid
        else:
            hi = mid

        if (hi - lo) < tol:
            return 0.5 * (lo + hi)

    # if not converged
    return 0.5 * (lo + hi)


def compute_radii_from_energy(E_blast_j: float, thresholds_psi=(2.0, 5.0, 10.0, 20.0)):
    """
    Top-level convenience function:
    - E_blast_j: energy coupled to blast (J)
    - returns dict: {psi: radius_m}
    """
    W = energy_to_kgTNT(E_blast_j)
    radii = {}
    for psi in thresholds_psi:
        R = find_radius_for_overpressure(W, psi)
        radii[psi] = R
    return W, radii


#feed this function energy in joules and eta (depends on terrain type)
#returns a list of tuples (psi, radius in meters) and seismic magnitude
def damage_coefficients_radii(kinetic_energy_joules, eta):
    blast_frac=0.3
    Wkg, radii = compute_radii_from_energy(kinetic_energy_joules*blast_frac)
    radius_pressures = []
    for psi,r in radii.items():
        if r is None:
            radius_pressures.append((psi, None))
        else:
            radius_pressures.append((psi, r))
    
    #creates a list of tuples (radius in meters, fatality rate, injury rate)
    coefficients_per_radius = []
    for r in radius_pressures:
        if r[0] == 2:
            coefficients_per_radius.append((r[1], 0.001, 0.5))
        elif r[0] == 5:
            coefficients_per_radius.append((r[1], 0.07, 0.6))  
        elif r[0] == 10:
            coefficients_per_radius.append((r[1], 0.3, 0.5))
        elif r[0] == 20:
            coefficients_per_radius.append((r[1], 0.7, 0.3))
    #returns a list of tuples (radius in meters, fatality rate, injury rate) and seismic magnitude
    return coefficients_per_radius,get_seismic_equivalent(kinetic_energy_joules, eta)

#CHEATSHEET:
# psi to kPa: multiply by 6.89476
# 2psi = light structural damage, windows shatter, lots of injuries but fatality rate 0.1%
# 5psi = moderate/heavy damage, wooden houses collapse, fatality rate ~1-10%
# 10psi = severe destruction, reinforced concrete falls, houses gone. Fatality rate ~10-50%
# 20psi = near total destruction, most buildings destroyed. Fatality rate ~50-90%




