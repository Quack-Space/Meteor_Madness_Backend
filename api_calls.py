import requests
import json

url = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"

def initiate_NEO():
    basic_params = "spkid,full_name,pdes,name,neo,pha,moid,e,a,q,i,om,w,ma,tp,per,n,ad,H,diameter,extent,GM,density,rot_per,albedo"
    display_params = ",source,producer,first_obs,last_obs,spec_T,spec_B"
    advanced_params = ",epoch,equinox,soln_date"

    keys = [
        "spkid", "full_name", "pdes", "name", "neo", "pha", "moid", "e", "a", "q", "i", 
        "om", "w", "ma", "tp", "per", "n", "ad", "H", "diameter", "extent", "GM", "density", 
        "rot_per", "albedo", "source", "producer", "first_obs", "last_obs", "spec_T", 
        "spec_B", "epoch", "equinox", "soln_date"
    ]

    params = {
        "fields": basic_params+display_params+advanced_params,
        "sb-group": "pha",
        "sb-kind": "a",
    }

    resp = requests.get(url, params=params)
    resp.raise_for_status()
    recived = resp.json()
    #print(recived)
    elem = recived["data"]
    dtype = lambda x: str(x).strip() if  not str(x).replace(".", "", 1).replace("-", "", 1).isdigit() else float(x)
    data = {l[0]: {keys[i]:  dtype(l[i]) for i in range(len(keys))} for l in elem}
    return data

def Neo_index():
    params = {"fields": "spkid","sb-group": "pha","sb-kind": "a",}
    resp = requests.get(url, params=params)
    data = resp.json()
    #data = json.load(data)
    index = [x[0] for x in data["data"]]
    print(index)
    print(len(index))
    return index


def NEO_catalog(db, start, offset):
    return db[start:start+offset]


def NEO_by_id(db, obj_id):
    asteroid = next((item for item in db if item["spkid"] == obj_id), None)
    return asteroid

