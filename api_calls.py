import requests
import json

## NEO db ##

def initiate_NEO():

    url = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"

    keys = [
        "spkid", "full_name", "pdes", "name", "neo", "pha", "moid", "e", "a", "q", "i", 
        "om", "w", "ma", "tp", "per", "n", "ad", "H", "diameter", "extent", "GM", "density", 
        "rot_per", "albedo", "source", "producer", "first_obs", "last_obs", "spec_T", 
        "spec_B", "epoch", "equinox", "soln_date"
    ]

    params = {
        "fields": ",".join(keys),
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
    url = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"
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


#################################################################################

#This code gets coordinates (lat,lon) and from those coordinates returns
#k, constant which if you multiply by the kinetic energy of the meteor at impact
#gives you the energy of the impact which u can convert to magnitude

#get elevation using etopo1
def get_elevation(lat, lon):
    url = f"https://api.opentopodata.org/v1/etopo1?locations={lat},{lon}"
    resp = requests.get(url)
    resp.raise_for_status()  # raises error if request failed
    data = resp.json()
    elev = data["results"][0]["elevation"]  # in meters
    return elev

#get rock type 👍
#if u touch this u dead i have no clue how this works
def get_rock_type(lat, lon): #officialy, it should be called "GLiM Class" not rock type
    url = "https://services8.arcgis.com/4KhTMTZ1x0f76DSg/arcgis/rest/services/GLiM_Niveau_I/FeatureServer/1/query"
    params = {
        "f": "json",
        "geometry": f"{lon},{lat}",          
        "geometryType": "esriGeometryPoint",
        "inSR": 4326,
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": "xx_Description,Litho,IDENTITY_",
        "returnGeometry": "false"
    }
    r = requests.get(url, params=params, timeout=20)
    r.raise_for_status()
    js = r.json()
    feats = js.get("features", [])
    if not feats:
        return None
    attrs = feats[0]["attributes"]
    # main label: "xx_Description" (e.g., "Unconsolidated sediments")
    return attrs.get("xx_Description") or attrs.get("Litho")


#costante di conversione energia cinetica e suoi possibili valori
k=0
K_LOOKUP = {
    "Unconsolidated sediments": 3e-4,
    "Siliciclastic sedimentary rocks": 5e-4,
    "Mixed sedimentary rocks": 5e-4,
    "Carbonate sedimentary rocks": 7e-4,
    "Metamorphics": 1e-3,
    "Acid/Intermediate/Basic plutonic rocks": 1e-3,
    "Acid/Intermediate/Basic volcanic rocks": 1e-3,
    "Pyroclastics": 7e-4,
    "Evaporites": 5e-4,
    "Ice and Glaciers": 5e-5,
    "Water Bodies": 1e-4,   
}

#ottieni k per determinate lat e lon
def get_k_constant(lat,lon):
    #altitudine punto
    elevation = get_elevation(lat, lon)
    # rock type 👍
    rock_type = get_rock_type(lat,lon)
    if elevation < 0:
        #it's water
        k=K_LOOKUP.get("Water Bodies")
    elif elevation > 0:
        #it's not water (unless it's a lake but genuinely. like there's no chance. let's be real)
        #tiziano says "what if it's a big fucking lake"
        #dev reply: nigga sybau
        k = K_LOOKUP.get(rock_type, 1e-3) 
    return k
