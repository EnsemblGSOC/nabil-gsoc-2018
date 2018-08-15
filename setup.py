"""
    This python script downloads the data and 
    installs the necessary python modules
"""


import requests
import pip
import os



data_url = "https://drive.google.com/uc?export=download&id=1X51bFX-WVbsp346byt36VW3kCUb6mPtb"
session = requests.Session()
response = session.get(data_url, params = { 'id' : id }, stream = True)

CHUNK_SIZE = 32768
with open(os.path.join('data_acquistion','data','db.db'), "wb") as f:
    for chunk in response.iter_content(CHUNK_SIZE):
        if chunk: 
            f.write(chunk)


pip.main(['install', 'Biopython'])
pip.main(['install', 'tqdm'])
pip.main(['install', 'Flask'])

