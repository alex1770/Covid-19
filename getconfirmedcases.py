from typing import Iterable, Dict, Union, List
from json import dumps
import sys,requests,json
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
    query_filters = [
        f"areaType=region",
        f"areaName=London"
    ]

    query_structure = {
        "date": "date",
        "name": "areaName",
        "code": "areaCode",
        "daily": "newCasesByPublishDate",
    }

    json_data = get_paginated_dataset(query_filters, query_structure)

    try:
        with open('confirmed.json','r') as fp:
            a=json.load(fp)
            olddate=a[0]['date']
    except:
        olddate='0000-00-00'

    if json_data[0]['date']<=olddate: sys.exit(1)
    
    with open('confirmed.json','w') as fp:
        json.dump(json_data,fp,indent=2)
    
    with open('confirmed.csv','w') as fp:
        pop=8831917# Use Zoe figure for London population
        start=0
        print("Date,London",file=fp)
        for x in json_data[::-1]:
            if x['daily']: start=1
            if start: print(x['date']+",","%7.1f"%(x['daily']/pop*1e6),file=fp)
