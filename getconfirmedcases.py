from typing import Iterable, Dict, Union, List
from json import dumps
import sys,requests,json,csv
from http import HTTPStatus

StructureType = Dict[str, Union[dict, str]]
FiltersType = Iterable[str]
APIResponseType = Union[List[StructureType], str]

def get_paginated_dataset(filters: FiltersType, structure: StructureType) -> APIResponseType:
    """
    Extracts paginated data by requesting all of the pages and combining the results.

    Parameters
    ----------
    filters: Iterable[str]
        API filters. See the API documentations for additional information.

    structure: Dict[str, Union[dict, str]]
        Structure parameter. See the API documentations for additional information.

    Returns
    -------
    Union[List[StructureType], str]
        Comprehensive list of dictionaries containing all the data for the given ``filters`` and ``structure``.
    """
    endpoint = "https://api.coronavirus.data.gov.uk/v1/data"

    api_params = {
        "filters": str.join(";", filters),
        "structure": dumps(structure, separators=(",", ":")),
        "format": "json"
    }

    data = list()

    page_number = 1

    while True:
        # Adding page number to query params
        api_params["page"] = page_number

        #print("Requesting page",page_number)
        response = requests.get(endpoint, params=api_params, timeout=10)

        if response.status_code >= HTTPStatus.BAD_REQUEST:
            raise RuntimeError(f'Request failed: {response.text}')
        elif response.status_code == HTTPStatus.NO_CONTENT:
            break

        current_data = response.json()
        page_data: List[StructureType] = current_data['data']
        
        data.extend(page_data)

        # The "next" attribute in "pagination" will be `None` when we reach the end.
        if current_data["pagination"]["next"] is None:
            break

        page_number += 1

    return data


if __name__ == "__main__":
    
    locs=[
        ("nation","England"),
        ("nation","Northern Ireland"),
        ("nation","Scotland"),
        ("nation","Wales"),
        ("region","East Midlands"),
        ("region","East of England"),
        ("region","London"),
        ("region","North East"),
        ("region","North West"),
        ("region","South East"),
        ("region","South West"),
        ("region","West Midlands"),
        ("region","Yorkshire and The Humber"),
    ]
    locs1=[y for (x,y) in locs]
    
    try:
        with open('confirmed.csv','r') as fp:
            reader=csv.reader(fp)
            for row in reader:
                pass
            olddate=row[0]
    except:
        olddate='0000-00-00'
        
    query_structure = {
        "date": "date",
        "name": "areaName",
        "code": "areaCode",
        "daily": "newCasesByPublishDate",
    }

    output={}
    dates=set()
    for (atype,aname) in locs:
      query_filters = [
        f"areaType="+atype,
        f"areaName="+aname
      ]
  
      json_data = get_paginated_dataset(query_filters, query_structure)
      if json_data[0]['date']<=olddate: sys.exit(1)
      #with open('confirmed.json','w') as fp: json.dump(json_data,fp,indent=2)

      output[aname]={}
      start=0
      for x in json_data[::-1]:
        if x['daily']: start=1
        if start: 
          dates.add(x['date'])
          output[aname][x['date']]=x['daily']
      
    with open('confirmed.csv','w') as fp:
      writer=csv.writer(fp)
      l=sorted(list(dates))
      writer.writerow(['Date']+locs1)
      #for loc in locs1: print(output[loc])
      for date in l:
        writer.writerow([date]+[output[loc].get(date,0) for loc in locs1])
