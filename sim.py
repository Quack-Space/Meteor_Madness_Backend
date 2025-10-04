import mehcanincs as m

def static_orbit(asteroid: dict):
    """Generate static orbit points from asteroid dict.
    Expects keys: 'a','e','i','om','w'. If any missing, returns empty list.
    """
    required = ["a","e","i","om","w"]
    if not asteroid or any(k not in asteroid for k in required):
        return []
    pts = m.generate_ellipse_points_shrf(
        asteroid["a"], asteroid["e"], i_deg=asteroid["i"],
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
    """Async placeholder full simulation.

    For now it just returns the static orbit data packaged as SimulationData after a tiny await.
    """
    import asyncio
    await asyncio.sleep(0.1)
    pts = static_orbit(data)
    traj = [{"position": {"x": p[0], "y": p[1], "z": p[2]}} for p in pts]
    return {
        "asteroid_id": str(data.get("spkid")),
        "asteroid_trajectory": traj,
        "earth_trajectory": [],
        "impact_estimate": {"will_impact": False},
    }