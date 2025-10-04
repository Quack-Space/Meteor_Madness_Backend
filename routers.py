from fastapi.middleware.cors import CORSMiddleware
from fastapi import HTTPException
from pydantic import BaseModel
from datetime import datetime, timedelta
from fastapi import FastAPI
from fastapi import Depends
from fastapi import status
import random
import asyncio
import uvicorn

#custom lib imports
import api_calls as api
import sim

#startup

app = FastAPI(
    title="QuackSpaceMeteorMadness",
    version="0.1"
)

status = 1;
r_count = 0;
startup_time = datetime.now()

async_results = {}


###START DB###
NEO = api.initiate_NEO()

####GET REQUESTS####

@app.get("/")
def root():
    uptime = timedelta(datetime.now(), startup_time)
    return {"status": status, "r_count": r_count, "uptime": uptime}

@app.get("/asteroids/")
def asteroids_catalog(start: int, offset: int):
    catalog = api.NEO_catalog(NEO, start, offset)
    return  {"asteroids": catalog}

@app.get("/asteroids/id/{obj_id}")
def asteroid_data(obj_id: str):
    data = api.NEO_by_id(obj_id)
    return {"data": data}

@app.get("/simulate/{obj_id}")
async def full_sim(obj_id: str, preview: bool):
    if preview:
        data = api.NEO_by_id(obj_id, job_id=None)
        result = sim.static_orbit(data)
        return result
    else:
        job_id = random.randint(00000, 99999)
        asyncio.create_task(sim.full_sim(data, job_id))
        return {"job_id": job_id}

@app.get("simulate/results/{job_id}")
def get_results(job_id: int):
    result = async_results[job_id]["data"]
    status = async_results[job_id]["status"]
    if status == True:
        del async_results[job_id]
    return {"status": status, "result": result}

####POST REQUESTS####

@app.post("/simulate_impact")
def impact_sim(data: dict):
    result = sim.impact(data)
    return result
    
uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)