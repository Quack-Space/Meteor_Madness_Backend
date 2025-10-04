import mehcanincs as m

def static_orbit(asteroid):
     pts = m.generate_ellipse_points_shrf(asteroid.a, asteroid.e, i_deg=asteroid.i,
                                          raan_deg=asteroid.om, argp_deg=asteroid.w,
                                          num_points=1000)
     return pts

def impact():
    pass

async def full_sim():
    pass