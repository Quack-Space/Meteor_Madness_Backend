import mehcanincs as m
import asyncio


def static_orbit(asteroid: dict):
    pts = m.generate_ellipse_points_shrf(asteroid["a"], asteroid["e"], i_deg=asteroid["i"],
                                        raan_deg=asteroid["om"], argp_deg=asteroid["w"],
                                         num_points=1000
                                         )
    return pts

def impact(data):
    #m, 2r, v, rho, alpha, latlon, terrain
    t_type = None
    diameter, depth = m.crater_dimensions_advanced(data["m"], data["v"], data["d"],
                                                   data["rho"], data["alpha"],
                                                   t_type)
    casualties = None
    sismic_magnitude = None

    return {"2r": diameter, "depth": depth}

async def full_sim(data: dict):
    pass