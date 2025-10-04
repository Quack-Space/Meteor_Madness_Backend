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

    def to_number(val):
        s = str(val).strip()
        # Accept plain int or float pattern but skip if leading zeros weird or too large to matter
        try:
            if s.replace('.', '', 1).replace('-', '', 1).isdigit():
                return float(s)
        except Exception:
            pass
        return None

    data = {}
    for row in elem:
        raw_id = str(row[0]).strip()
        # Ensure no trailing .0 artifacts (e.g., '2000433.0')
        if raw_id.endswith('.0') and raw_id.replace('.0', '').isdigit():
            raw_id = raw_id[:-2]
        entry = {}
        for i, k in enumerate(keys):
            val = row[i]
            if k == 'spkid':
                entry[k] = raw_id  # keep canonical id string
                continue
            num_val = to_number(val)
            entry[k] = num_val if num_val is not None else (str(val).strip() if val not in (None, '') else None)
        data[raw_id] = entry
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


def NEO_catalog(db, start: int, offset: int):
    """Return a slice of the NEO catalog.

    db: dict keyed by spkid -> asteroid data dict
    start: starting index in the ordered key list
    offset: number of records to return
    """
    if offset <= 0:
        return []
    keys = list(db.keys())
    slice_keys = keys[start:start+offset]
    return [db[k] for k in slice_keys]

def NEO_by_id(db, obj_id):
    """Lookup asteroid by spkid directly in dict (returns None if missing)."""
    return db.get(obj_id)

