from datetime import datetime
from fastapi import FastAPI, HTTPException, Response, status
from fastapi.middleware.cors import CORSMiddleware
import random
import asyncio
from typing import Dict, Any, List

#custom lib imports
import api_calls as api
import sim

#startup

app = FastAPI(
    title="QuackSpaceMeteorMadness",
    version="0.1"
)

service_status = 1
request_count = 0
startup_time = datetime.utcnow()

async_results = {}


###START DB###
NEO = api.initiate_NEO()

# CORS (frontend likely on Vite default http://localhost:5173)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

####GET REQUESTS####

@app.get("/")
def root():
    uptime = (datetime.utcnow() - startup_time).total_seconds()
    return {"status": service_status, "request_count": request_count, "uptime_seconds": uptime}

def _map_asteroid(raw: Dict[str, Any]) -> Dict[str, Any]:
    # Provide safe numeric parsing helpers
    def num(val):
        return val if isinstance(val, (int, float)) else None
    diameter = num(raw.get("diameter"))  # Already in km per NASA small-body db (check upstream)
    density = num(raw.get("density"))  # Might be None often
    # Derive mass if both diameter & density exist
    mass = None
    if diameter and density and diameter > 0 and density > 0:
        radius_m = (diameter * 1000) / 2.0
        volume_m3 = 4/3 * 3.141592653589793 * radius_m**3
        mass = volume_m3 * density  # kg
    rotation_period = num(raw.get("rot_per"))  # hours already? Assume hours
    mapped = {
        "id": str(raw.get("spkid")),
        "name": raw.get("name") or raw.get("full_name") or str(raw.get("spkid")),
        "diameter": diameter or 0,
        "density": density or 0,
        "velocity": None,  # Not provided by current dataset
        "mass": mass,
        "albedo": num(raw.get("albedo")),
        "rotation_period": rotation_period,
        "spectral_type": raw.get("spec_T") or raw.get("spec_B"),
    }
    return mapped

@app.get("/asteroids", response_model=None)
def asteroids_catalog(start: int = 0, limit: int = 200):
    if limit <= 0:
        raise HTTPException(status_code=400, detail="limit must be > 0")
    keys: List[str] = list(NEO.keys())
    slice_keys = keys[start:start+limit]
    items = [_map_asteroid(NEO[k]) for k in slice_keys]
    return {
        "asteroids": items,
        "total": len(keys),
        "limit": limit,
        "offset": start,
    }

@app.get("/asteroids/{obj_id}")
def asteroid_data(obj_id: str):
    raw = api.NEO_by_id(NEO, obj_id)
    if raw is None:
        raise HTTPException(status_code=404, detail="Asteroid not found")
    return _map_asteroid(raw)

@app.post("/simulate/{obj_id}")
async def simulate(obj_id: str, preview: bool = True):
    raw = api.NEO_by_id(NEO, obj_id)
    if raw is None:
        raise HTTPException(status_code=404, detail="Asteroid not found")
    if preview:
        points = sim.static_orbit(raw)
        # Map points list of (x,y,z) => trajectory object list
        traj = [{"position": {"x": p[0], "y": p[1], "z": p[2]}} for p in points]
        return {
            "asteroid_id": obj_id,
            "asteroid_trajectory": traj,
            "earth_trajectory": [],
            "impact_estimate": {"will_impact": False},
        }

    job_id = random.randint(10000, 99999)
    async_results[job_id] = {"done": False, "data": None}

    async def run_full():
        try:
            res = await sim.full_sim(raw)
            # Expect res already shaped; if None produce placeholder
            if res is None:
                points = sim.static_orbit(raw)
                res = {
                    "asteroid_id": obj_id,
                    "asteroid_trajectory": [{"position": {"x": p[0], "y": p[1], "z": p[2]}} for p in points],
                    "earth_trajectory": [],
                    "impact_estimate": {"will_impact": False},
                }
            async_results[job_id] = {"done": True, "data": res}
        except Exception as e:
            async_results[job_id] = {"done": True, "data": {"error": str(e)}}

    asyncio.create_task(run_full())
    return {"job_id": str(job_id)}

@app.get("/reports/{job_id}")
def get_report(job_id: str, response: Response):
    # job ids stored as int keys originally -> unify string storage
    # Ensure key checks for both str/int
    job = async_results.get(job_id)
    if job is None and job_id.isdigit():
        job = async_results.get(int(job_id))
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if not job.get("done"):
        response.status_code = status.HTTP_202_ACCEPTED
        return {"status": "processing"}
    data = job.get("data")
    # Remove once delivered
    try:
        del async_results[job_id]
    except KeyError:
        try:
            del async_results[int(job_id)]
        except Exception:
            pass
    return data

####POST REQUESTS####

@app.post("/simulate_impact")
def impact_sim(data: dict):
    result = sim.impact(data)
    return result
    
if __name__ == "__main__":
    import uvicorn
    # Bind explicitly to localhost for dev convenience
    uvicorn.run("routers:app", host="127.0.0.1", port=8000, reload=True)
