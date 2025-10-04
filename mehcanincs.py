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
