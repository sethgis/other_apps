import ee
import datetime
import folium
import params
import os
import subprocess

# Import the matplotlib.pyplot module.
import matplotlib.pyplot as plt
# import matplotlib.pyplot as plt
import numpy as np
import requests
# import pillow
from PIL import Image
import shutil
import urllib
import numpy as np
import matplotlib.pyplot as plt
import plotly_express as px

import matplotlib.pyplot as plt
import numpy as np
!pip install pandas
import pandas as pd
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
import plotly.io as pio

try:
    import geemap
except ImportError:
    print('geemap package not installed. Installing ...')
    subprocess.check_call(["python", '-m', 'pip', 'install', 'geemap'])


MODEL = os.path.dirname(os.path.abspath("__file__"))

ppk_file = os.path.join(MODEL, "service_account.json")

service_account = "seth-311@seth-.iam.gserviceaccount.com"

private_key = ppk_file

credentials = ee.ServiceAccountCredentials(service_account, private_key)
ee.Initialize(credentials)


""" ACCESSING AREA OF INTEREST"""

# table = ee.FeatureCollection('projects/ee-snyawacha/assets/CDR_BOUNDARIES/KE_Admin_Level_1')
# table1 = table.geometry()


AOI = ee.FeatureCollection("projects/mlearnings/assets/KITUI_UPDATE_2024").filter(ee.Filter.inList('AEZ', ['6'])).geometry()

title = 'AEZ_6'
# table2 = table.filter(ee.Filter.inList('CLUSTER_1', ['7'])).geometry()







"""
COMPUTING THE LONG TERM AVERAGE RAINFALL FOR A GIVEN LOCATION
 """

months = ee.List.sequence(250, 365, 5)

sst = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')\
    .filterDate(ee.Date('2018-01-01'), ee.Date('2023-12-30'))

""" MAKE A FUNCTION TO COMPUTE A DAILY AVERAGE RAINFALL OVER THE PAST FIVE YEARS"""

def function2(m):
    return sst.filter(ee.Filter.calendarRange(m,m, 'day_of_year'))\
    .mean().set('day_of_year', m).set('system:time_start', sst.get('system:time_start'))
    
""" CALLING THE FUNCTION AND PRINTING OUT THE OUTPUTS"""
LTA = ee.ImageCollection.fromImages(months.map(function2))

LTA


""" ACCESSING THE DAILY RAINFALL PATTERN FOR THE YEAR 2023"""

DAILY = ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD')\
        .filterDate(ee.Date('2023-09-01'), ee.Date('2023-12-30'))
 
    

BANDS = DAILY.toBands()
LTA_BANDS = LTA.toBands().set('system:time_start', sst.get('system:time_start'))
# AOI = ee.FeatureCollection('projects/ee-snyawacha/assets/Machakos')
LTA

# BANDS
# LTA_BANDS





series = BANDS.reduceRegions(collection=AOI,
                                      reducer=ee.Reducer.mean(),
                                      scale=100)

dic_series = series.getInfo()


dic_series

EVI_anom = np.array(list(dic_series['features'][0]['properties'].values()))
# print(EVI_anom)
array1 = pd.DataFrame(EVI_anom).reset_index()




series2 = LTA_BANDS.reduceRegions(collection=AOI,
                                      reducer=ee.Reducer.mean(),
                                      scale=100)

dic_series2 = series2.getInfo()


dic_series2
EVI_anom2 = np.array(list(dic_series2['features'][0]['properties'].values()))
# print(EVI_anom)
array2 = pd.DataFrame(EVI_anom2).reset_index()

array2

array1


import datetime
import pandas as pd
 
# initializing date
test_date = datetime.datetime.strptime("01-09-2023", "%d-%m-%Y")
 

K = 120
 
date_generated = pd.date_range(test_date, periods=K)



date_gen = date_generated[::5]

date_frame = (pd.DataFrame(date_gen)).reset_index()
date_frame

merged = pd.merge(pd.merge(array1,array2,on='index'),date_frame,on='index')
merged




import plotly_express as px

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
  
merged2 = merged
# merged2.columns()
merged3 = merged2.rename(columns = {'0_x':'Daily', '0_y':'LTA',0:'Date' })
merged3

merged3

df2 = merged3.iloc[:-2 , :]
df2

####Drawing the graph of rainfall plot for the specific AEZ

# df2['Daily']=df2(['Daily']).astype(float)
# df2['LTA']=df2(['LTA']).astype(float)


# df.drop(df[df['Fee'] >= 24000].index, inplace = True)


df = df2

# df2['LTA'].plot(kind='line', marker='*', color='black', ms=10)


from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import plotly


# import plotly.io as pio
# pio.renderers.default = "sphinx_gallery"

# data

# set up plotly figure
fig = make_subplots(1,2)

# add first bar trace at row = 1, col = 1
fig.add_trace(go.Bar(x=df['Date'], y=df['Daily'],
                     name='Daily Rainfall Amount',
                     marker_color = 'blue',
                     marker_line_color='rgb(8,48,107)',
                     marker_line_width=2),
              row = 1, col = 1)


fig.add_trace(go.Line(x=df['Date'], y=df['LTA'],
                     name='Long Term Average Rainfall',
                     marker_color = 'red',
                     marker_line_color='rgb(8,48,107)',
                    marker_line_width=3),
              row = 1, col = 1)

# fig.show()

# fig.update_layout(
#     autosize=False,
#     width=2000,
#     height=800,)
fig.update_layout(
    title= title,
    legend=dict(x=0, y=1, traceorder='normal', orientation='h'),
    autosize=False,
    width=2000,
    height=800,
)

# fig.show()

# pio.write_image(fig, file=title, format="png")

# print("Image Succeffully Saved")
# print(fig)


