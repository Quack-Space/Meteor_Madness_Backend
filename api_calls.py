import requests

url = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"
page=1
limit=5000
all_rows=[]

params = {
    "fields": "full_name,spkid,neo,pha,moid,H",
    "limit": limit,
    "sb-group": "pha",
}

resp = requests.get(url, params=params)
resp.raise_for_status()
data = resp.json()
