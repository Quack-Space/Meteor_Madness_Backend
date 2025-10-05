import mehcanincs as m
import api_calls as api
from typing import Dict, Any, Tuple, List
from datetime import datetime, timezone
import math
import time

AU_METERS = 149_597_870_700.0  # 1 AU in meters

# ---------------------------------------------------------------------------
# Physics helper utilities (intentionally lightweight: NO n-body / perturbations)
# ---------------------------------------------------------------------------

def period_and_mean_motion(a_m: float, mu: float) -> Tuple[float, float]:
    """Return (orbital_period_seconds, mean_motion_rad_s)."""
    T = 2.0 * math.pi * math.sqrt(a_m ** 3 / mu)
    n = math.sqrt(mu / a_m ** 3)
    return T, n

def normalize_angle_rad(theta: float) -> float:
    twopi = 2.0 * math.pi
    return theta % twopi

def anomalies_from_position_sequence(positions: List[List[float]], e: float) -> Tuple[List[float], List[float], List[float]]:
    """Derive true, eccentric, and mean anomalies for a planar ellipse sample.

    Assumes orbit lies in (approximately) a single plane aligned with XY (current simplified case).
    This avoids modifying mechanics code to return nu directly. Adequate for hackathon scope.
    Returns lists of (nu_rad, E_rad, M_rad).
    """
    if e < 0.0 or e >= 1.0:
        return [], [], []
    true_list: List[float] = []
    ecc_list: List[float] = []
    mean_list: List[float] = []
    sqrt_fac = math.sqrt((1 - e) / (1 + e)) if e < 1 else 0.0
    for p in positions:
        x, y = p[0], p[1]
        nu = math.atan2(y, x)
        # E from nu
        if abs(e) < 1e-12:
            E = nu
        else:
            tan_half_nu = math.tan(nu / 2.0)
            # Guard huge tan values; fallback to nu if overflow
            try:
                tan_half_E = sqrt_fac * tan_half_nu
                E = 2.0 * math.atan(tan_half_E)
            except Exception:
                E = nu
            E = normalize_angle_rad(E)
        M = E - e * math.sin(E)
        M = normalize_angle_rad(M)
        true_list.append(nu)
        ecc_list.append(E)
        mean_list.append(M)
    return true_list, ecc_list, mean_list

def compute_closest_approach(p1: List[List[float]], p2: List[List[float]]) -> Tuple[int, float]:
    """Return (index, distance_m) of minimum separation for synchronized samples."""
    if not p1 or not p2 or len(p1) != len(p2):
        return -1, float('nan')
    
    try:
        d = lambda a, b: (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2
        imin, dmin = min([(i, d(a,b)) for i, (a,b) in enumerate(zip(p1, p2))], key=lambda x: x[1])
    except ValueError:
            dmin = float('nan')
            imin = -1
    return imin, math.sqrt(dmin) if imin >= 0 else float('nan')


def earth_spin_metadata() -> dict:
    """Return Earth axial tilt + rotation model constants."""
    return {
        "obliquity_deg": 23.439281,          # mean obliquity (approx, J2000)
        "rotation_period_s": 86164.0905,     # sidereal day
        "prime_meridian_rad_at_epoch": 0.0   # reference; front-end can advance linearly
    }


def static_orbit(asteroid: dict):
    # Catalog gives 'a' in AU. Convert to meters for mechanics routine which expects meters.
    try:
        a_au = float(asteroid["a"])
    except (KeyError, TypeError, ValueError):
        return []
    e = float(asteroid.get("e", 0.0) or 0.0)
    i_deg = float(asteroid.get("i", 0.0) or 0.0)
    raan_deg = float(asteroid.get("om", 0.0) or 0.0)
    argp_deg = float(asteroid.get("w", 0.0) or 0.0)
    a_m = a_au * AU_METERS
    pts = m.generate_ellipse_points_shrf(a_m, e, i_deg=i_deg,
                                         raan_deg=raan_deg, argp_deg=argp_deg,
                                         num_points=1000)
    # Provide simple orbital metadata for diagnostics (AU based)
    q_au = a_au * (1 - e)
    Q_au = a_au * (1 + e)
    return {
        "points": pts,
        "orbit_meta": {
            "a_au": a_au,
            "e": e,
            "q_au": q_au,
            "Q_au": Q_au,
            "i_deg": i_deg,
            "raan_deg": raan_deg,
            "argp_deg": argp_deg,
        }
    }

async def full_sim(data: dict):
    """Generate raw arrays for asteroid full simulation (no per-point dicts).

    Returns a JSON-serializable dict with ONLY array primitives so the frontend can
    transform as it likes:
    {
       "asteroid_id": str,
       "asteroid_positions": [[x,y,z], ...],
       "asteroid_velocities": [[vx,vy,vz], ...],
       "timestamps": [t0, t1, ...],              # seconds since periapsis
       "earth_position": [x,y,z],                # current heliocentric Earth position (approx)
       "orbit_meta": { a_au, e, q_au, Q_au, i_deg, raan_deg, argp_deg, period_seconds }
    }
    """
    asteroid: Dict[str, Any] = data
    try:
        a_au = float(asteroid.get("a"))
    except (TypeError, ValueError):
        return {"error": "missing_or_invalid_a"}
    e = float(asteroid.get("e", 0.0) or 0.0)
    i_deg = float(asteroid.get("i", 0.0) or 0.0)
    raan_deg = float(asteroid.get("om", 0.0) or 0.0)
    argp_deg = float(asteroid.get("w", 0.0) or 0.0)

    a_m = a_au * AU_METERS
    positions, velocities, times = m.generate_ellipse_timed_points_shrf(
        a_m, e, i_deg=i_deg, raan_deg=raan_deg, argp_deg=argp_deg, num_points=1000
    )

    # Flatten rounding already applied inside mechanics.
    asteroid_positions = [[p[0], p[1], p[2]] for p in positions]
    asteroid_velocities = [[v[0], v[1], v[2]] for v in velocities]
    timestamps = times  # list of floats

    # Orbital metadata (AU based) + period (seconds)
    q_au = a_au * (1 - e)
    Q_au = a_au * (1 + e)
    period_seconds, mean_motion = period_and_mean_motion(a_m, m.mu_Sun)

    # Generate Earth trajectory sampled at same count using SAME index alignment
    earth_positions, earth_velocities, earth_times = _sample_earth_orbit(len(times))
    earth_a_m = AU_METERS
    earth_e = 0.0167
    earth_period_seconds, earth_mean_motion = period_and_mean_motion(earth_a_m, m.mu_Sun)

    # Compute anomalies (simplified planar assumption)
    asteroid_true, asteroid_ecc, asteroid_mean = anomalies_from_position_sequence(asteroid_positions, e)
    earth_true, earth_ecc, earth_mean = anomalies_from_position_sequence(earth_positions, earth_e)

    # Shared normalized progress 0..1 for interpolation convenience
    n_samples = len(timestamps)
    progress = [i / (n_samples - 1) if n_samples > 1 else 0.0 for i in range(n_samples)]

    # Closest approach scan (synchronized samples)
    ca_index, ca_distance = compute_closest_approach(asteroid_positions, earth_positions)

    # Provide epoch reference (POSIX seconds "now") so frontend can map absolute time if desired
    epoch_now = time.time()
    # Absolute timestamps for asteroid = epoch_now + (timestamps since perihelion interpreted as synthetic)
    asteroid_absolute = [epoch_now + t for t in timestamps]
    earth_absolute = [epoch_now + (earth_period_seconds * p) for p in progress]

    # Response (backwards compatible keys preserved + enriched metadata)
    return {
        "epoch": epoch_now,
        "units": {"length": "m", "time": "s"},
        "mu_sun": m.mu_Sun,
        "asteroid_id": str(asteroid.get("spkid") or asteroid.get("id") or "unknown"),
        # Original fields (maintain):
        "asteroid_positions": asteroid_positions,
        "asteroid_velocities": asteroid_velocities,
        "timestamps": timestamps,  # original non-uniform (seconds since periapsis)
        "earth_positions": earth_positions,
        # New / enriched fields:
        "earth_velocities": earth_velocities,
        "progress": progress,  # unified normalized timeline
        "asteroid_absolute_timestamps": asteroid_absolute,
        "earth_absolute_timestamps": earth_absolute,
        "asteroid_orbit": {
            "a_au": a_au,
            "e": e,
            "q_au": q_au,
            "Q_au": Q_au,
            "i_deg": i_deg,
            "raan_deg": raan_deg,
            "argp_deg": argp_deg,
            "period_seconds": period_seconds,
            "mean_motion_rad_s": mean_motion,
            "anomalies": {
                "true_anomaly_rad": asteroid_true,
                "eccentric_anomaly_rad": asteroid_ecc,
                "mean_anomaly_rad": asteroid_mean,
            }
        },
        "earth_orbit": {
            "a_au": 1.0,
            "e": earth_e,
            "period_seconds": earth_period_seconds,
            "mean_motion_rad_s": earth_mean_motion,
            "anomalies": {
                "true_anomaly_rad": earth_true,
                "eccentric_anomaly_rad": earth_ecc,
                "mean_anomaly_rad": earth_mean,
            },
            "spin": earth_spin_metadata(),
        },
        "closest_approach": {
            "index": ca_index,
            "distance_m": ca_distance,
            "progress": progress[ca_index] if ca_index >= 0 else None,
            "asteroid_time_s": timestamps[ca_index] if ca_index >= 0 else None,
            "earth_time_s": earth_times[ca_index] if ca_index >= 0 else None,
        },
        # Legacy alias retained for frontend backward compatibility
        "orbit_meta": {
            "a_au": a_au,
            "e": e,
            "q_au": q_au,
            "Q_au": Q_au,
            "i_deg": i_deg,
            "raan_deg": raan_deg,
            "argp_deg": argp_deg,
            "period_seconds": period_seconds,
        }
    }


def _earth_current_position_heliocentric() -> Tuple[float, float, float]:
    """Approximate current Earth heliocentric ecliptic position (meters).

    Simplified: inclination, RAAN ~ 0; argument of perihelion ignored; solves Kepler's
    equation relative to an approximate perihelion epoch.
    """
    # Orbital elements (simplified)
    a_m = AU_METERS
    e = 0.0167
    # Reference perihelion (approx 2025 Jan 4 00:00 UTC)
    ref_peri = datetime(2025, 1, 4, 0, 0, tzinfo=timezone.utc)
    now = datetime.now(timezone.utc)
    dt = (now - ref_peri).total_seconds()
    # Orbital period (sidereal year)
    T = 365.256363004 * 86400.0
    M = 2.0 * math.pi * ((dt % T) / T)
    # Solve Kepler's equation for E
    E = M
    for _ in range(8):
        f = E - e * math.sin(E) - M
        fp = 1 - e * math.cos(E)
        E -= f / fp
    # True anomaly
    nu = 2.0 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2.0), math.sqrt(1 - e) * math.cos(E / 2.0))
    r = a_m * (1 - e * math.cos(E))
    x = r * math.cos(nu)
    y = r * math.sin(nu)
    z = 0.0
    return round(x, 1), round(y, 1), round(z, 1)


def _sample_earth_orbit(n: int) -> Tuple[list, list, list]:
    """Sample an approximate Earth orbit (positions, velocities, times) with n samples.

    Simplified: use same eccentric approximation as _earth_current_position_heliocentric
    but spread uniformly in mean anomaly over one sidereal year.
    """
    if n <= 0:
        return [], [], []
    a_m = AU_METERS
    e = 0.0167
    T = 365.256363004 * 86400.0
    positions = []
    velocities = []
    times = []
    mu = m.mu_Sun
    for idx in range(n):
        M = 2.0 * math.pi * (idx / n)
        # Newton solve for E
        E = M
        for _ in range(6):
            f = E - e * math.sin(E) - M
            fp = 1 - e * math.cos(E)
            E -= f / fp
        nu = 2.0 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2.0), math.sqrt(1 - e) * math.cos(E / 2.0))
        r = a_m * (1 - e * math.cos(E))
        x = r * math.cos(nu)
        y = r * math.sin(nu)
        z = 0.0
        # perifocal velocity magnitude components (simplified planar)
        p = a_m * (1 - e * e)
        h = math.sqrt(mu * p)
        vx_pf = -mu / h * math.sin(nu)
        vy_pf = mu / h * (e + math.cos(nu))
        # For planar orbit, perifocal == inertial XY
        positions.append([round(x, 1), round(y, 1), 0.0])
        velocities.append([round(vx_pf, 1), round(vy_pf, 1), 0.0])
        times.append(round(T * (idx / n), 1))
    return positions, velocities, times




def deaths_and_injured(kinetic_energy, crater_d, eta, lat, lon):
    death_count = 0
    injured_count = 0

    #crater diameter, everybody dead
    pop_crater = api.pop_within_radius_ghs(lat,lon,math.ceil(crater_d/2))
    pop_crater_normalized = pop_crater*(crater_d/math.ceil(crater_d)) #statistic
    death_count += pop_crater_normalized

    #get damage radius with associated fatality and injury rates
    psi_radius = m.damage_coefficients_radii(kinetic_energy, eta)

    #now with the radii
    pop_20psi = api.pop_within_radius_ghs(lat,lon,math.ceil(psi_radius[0][3][0]))-pop_crater_normalized
    dead_20psi = pop_20psi*psi_radius[0][3][1]
    injured_20psi = (pop_20psi-dead_20psi)*psi_radius[0][3][2]
    death_count += dead_20psi
    injured_count += injured_20psi

    pop10psi = api.pop_within_radius_ghs(lat,lon,math.ceil(psi_radius[0][2][0]))-pop_crater_normalized-pop_20psi
    dead_10psi = pop10psi*psi_radius[0][2][1]
    injured_10psi = (pop10psi-dead_10psi)*psi_radius[0][2][2]   
    death_count += dead_10psi
    injured_count += injured_10psi

    pop5psi = api.pop_within_radius_ghs(lat,lon,math.ceil(psi_radius[0][1][0]))-pop_crater_normalized-pop_20psi-pop10psi
    dead_5psi = pop5psi*psi_radius[0][1][1]
    injured_5psi = (pop5psi-dead_5psi)*psi_radius[0][1][2]
    death_count += dead_5psi
    injured_count += injured_5psi

    pop2psi = api.pop_within_radius_ghs(lat,lon,math.ceil(psi_radius[0][0][0]))-pop_crater_normalized-pop_20psi-pop10psi-pop5psi
    dead_2psi = pop2psi*psi_radius[0][0][1]
    injured_2psi = (pop2psi-dead_2psi)*psi_radius[0][0][2]
    death_count += dead_2psi
    injured_count += injured_2psi

    #list of radii for infogrpahics purposes
    radii = [crater_d/2, psi_radius[0][3][0], psi_radius[0][2][0], psi_radius[0][1][0], psi_radius[0][0][0]]

    #returns: [0], number of dead, [1] number of injured, [2] list of radii (for infogrpahics purposes)
    return math.ceil(death_count), math.ceil(injured_count), radii
    

def impact(data):
    #m, 2r, v, rho, alpha, latlon, terrain
    """
    0: crater
    1: pressure wave 20psi
    2: pressure wave 10psi
    3: pressure wave 5psi
    4: pressure wave 2psi
    """
    eta, t_type = api.get_terrain_characteristics(data["lat"], data["lon"])
    rho = data["m"]/data["v"]
    c_diameter, c_depth, Ek_impact, m_abl = m.crater_dimensions_advanced(data["m"], data["v"], data["d"],
                                                   rho, data["alpha"], t_type)
    sismic_magnitude = m.get_seismic_equivalent(Ek_impact, eta)
    deads, injured, radii = deaths_and_injured(Ek_impact, c_diameter, eta, data["lat"], data["lon"])
    

    return {"ablation": m_abl, "Ek_impact": Ek_impact, "magnitude": sismic_magnitude,
            "deaths": deads, "injured": injured, "radii": radii,
            "crater": {"diameter": c_diameter, "depth": c_depth}}
