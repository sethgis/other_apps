import os
import json
import geopandas as gpd
import rasterio
from rasterio.mask import mask
import numpy as np

# you must first reproject the geotiff to version that matches the shapefiles, preferebaly epsg:4326

def fetch_and_save_yield():
    # Load the main GeoJSON with many farms
    with open('/Users/sethnyawacha/Desktop/APIs/FARM_BOUNDARU_DRAWER/downloaded/JOINED_15K_LAYER.geojson') as f:
        data = json.load(f)

    years = ['2019', '2020', '2021', '2022', '2023']  # Correspond to the bands in each GeoTIFF
    
    # Paths to the multi-band GeoTIFFs for AGBD, yield, and SOC
    geotiffs = {
        'season_A': {
            'agbd': "/Users/sethnyawacha/Desktop/APIs/FARM_BOUNDARU_DRAWER/stacked_2/PROJ_AGBD_SEASON_A.tif",
            'yield': "/Users/sethnyawacha/Desktop/APIs/FARM_BOUNDARU_DRAWER/stacked_2/PROJ_YIELD_SEASON_A.tif",
            "spi": "/Users/sethnyawacha/Desktop/APIs/FARM_BOUNDARU_DRAWER/stacked_2/resampled_file_season_A.tif",
            'soc':"/Volumes/SanDisk/SAN_ASOC_2019_2023_SEASON_A_MERGED.tif"
            
        },
        'season_B': {
            'agbd': "/Users/sethnyawacha/Desktop/APIs/FARM_BOUNDARU_DRAWER/stacked_2/PROJ_AGBD_SEASON_B.tif",
            'yield': "/Users/sethnyawacha/Desktop/APIs/FARM_BOUNDARU_DRAWER/stacked_2/PROJ_YIELD_SEASON_B.tif",
            "spi": "/Users/sethnyawacha/Desktop/APIs/FARM_BOUNDARU_DRAWER/stacked_2/resampled_file_season_B.tif",
            'soc':"/Volumes/SanDisk/SAN_ASOC_2019_2023_SEASON_B_MERGED.tif"
        }
    }

    for index, feature in enumerate(data['features']):
        feature_id = feature['properties'].get('FID')
        individual_geojson = {
            "geojson": {
                "type": "FeatureCollection",
                "features": [feature]
            }
        }

        # Initialize storage for statistics values for each season
        values = {
            'season_A': {year: {'agbd': None, 'yield': None, 'spi': None, 'soc': None} for year in years},
            'season_B': {year: {'agbd': None, 'yield': None, 'spi': None, 'soc': None} for year in years}
        }

        # Process each season
        for season, paths in geotiffs.items():
            for geotiff_type, geotiff_path in paths.items():
                try:
                    # Open the multi-band GeoTIFF
                    with rasterio.open(geotiff_path) as src:
                        boundary_data = gpd.GeoDataFrame.from_features(individual_geojson["geojson"]["features"])
                        boundary_geom = boundary_data.geometry.values[0]
                        boundary_geojson = boundary_geom.__geo_interface__

                        # Mask the raster with the farm boundary
                        masked_image, _ = mask(src, [boundary_geojson], crop=True, nodata=np.nan)

                        # Loop through the bands (each band corresponds to a year)
                        for band_index in range(masked_image.shape[0]):
                            year = years[band_index]

                            # Extract valid pixel values for the current band (year)
                            pixel_values = masked_image[band_index][~np.isnan(masked_image[band_index])]
                            mean_value = float(np.mean(pixel_values)) if pixel_values.size > 0 else 0

                            # Store the mean value in the appropriate dictionary
                            values[season][year][geotiff_type] = mean_value

                except FileNotFoundError:
                    print(f"GeoTIFF file not found: {geotiff_path}. Skipping...")

        # Combine the statistics into the feature's properties
        combined_data = {
            'type': 'Feature',
            'geometry': feature['geometry'],
            'properties': {
                'feature_id': feature_id,
                'season_A': {
                    'agbd': {year: values['season_A'][year]['agbd'] for year in years},
                    'yield': {year: values['season_A'][year]['yield'] for year in years},
                    'spi': {year: values['season_A'][year]['spi'] for year in years},
                    'soc': {year: values['season_A'][year]['soc'] for year in years}
                },
                'season_B': {
                    'agbd': {year: values['season_B'][year]['agbd'] for year in years},
                    'yield': {year: values['season_B'][year]['yield'] for year in years},
                    'spi': {year: values['season_B'][year]['spi'] for year in years},
                    'soc': {year: values['season_B'][year]['soc'] for year in years}
                },
                'years': years
            }
        }

        # Save the results in a GeoJSON file for this feature
        geojson_filename = f"/Volumes/SanDisk/FARMS/combined_{index}.geojson"
        with open(geojson_filename, 'w') as geojson_file:
            json.dump(combined_data, geojson_file)

        print(f"Saved combined GeoJSON data for feature {index} in {geojson_filename}")

# Call the function
fetch_and_save_yield()
