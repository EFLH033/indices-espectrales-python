
"""
Created on Wed Jul 16 16:06:57 2025

@author: Familia
"""

import rasterio
import numpy as np
import matplotlib.pyplot as plt
import os
import geopandas as gpd
from rasterio.enums import Resampling
from rasterio.warp import reproject, calculate_default_transform
from matplotlib.colors import LinearSegmentedColormap
import rasterio.crs
from rasterio.mask import mask
from skimage import exposure

# --- Función para Cargar GeoJSON ---
def load_geojson(geojson_file_path):
    """
    Carga un archivo GeoJSON en un GeoDataFrame.
    """
    try:
        gdf = gpd.read_file(geojson_file_path)
        print(f"--- GeoJSON de Polígonos Cargado: {os.path.basename(geojson_file_path)} ---")
        print(f"Número de polígonos encontrados: {len(gdf)}")
        # Asegurarse de que el GeoJSON esté en WGS84
        if gdf.crs != 'EPSG:4326':
            gdf = gdf.to_crs('EPSG:4326')
            print(f"El GeoJSON se cargó exitosamente en EPSG:4326.")
        else:
            print(f"El GeoJSON se cargó exitosamente en {gdf.crs}.")
        return gdf
    except Exception as e:
        print(f"Error al cargar el GeoJSON: {e}")
        return None

# --- DEFINICIÓN DE FUNCIONES DE ÍNDICES ---
def calcular_ndvi(nir_band, red_band):
    """
    Calcula el Índice de Vegetación de Diferencia Normalizada (NDVI).
    """
    nir_band_float = nir_band.astype(np.float32)
    red_band_float = red_band.astype(np.float32)
    numerator = nir_band_float - red_band_float
    denominator = nir_band_float + red_band_float
    # Manejo de división por cero: donde el denominador es 0, el resultado es 0.
    ndvi = np.where(denominator == 0, 0, numerator / denominator)
    return ndvi

def calcular_ndmi(nir_band, swir_band):
    """
    Calcula el Índice de Humedad de Diferencia Normalizada (NDMI).
    """
    nir_band_float = nir_band.astype(np.float32)
    swir_band_float = swir_band.astype(np.float32)
    numerator = nir_band_float - swir_band_float
    denominator = nir_band_float + swir_band_float
    ndmi = np.where(denominator == 0, 0, numerator / denominator)
    return ndmi

def calcular_ndwi(green_band, nir_band):
    """
    Calcula el Índice de Diferencia Normalizada del Agua (NDWI).
    """
    green_band_float = green_band.astype(np.float32)
    nir_band_float = nir_band.astype(np.float32)
    numerator = green_band_float - nir_band_float
    denominator = green_band_float + nir_band_float
    ndwi = np.where(denominator == 0, 0, numerator / denominator)
    return ndwi

def calcular_evi(nir_band, red_band, blue_band):
    """
    Calcula el Índice de Vegetación Mejorado (EVI).
    """
    nir_band_float = nir_band.astype(np.float32)
    red_band_float = red_band.astype(np.float32)
    blue_band_float = blue_band.astype(np.float32)
    L = 1.0
    C1 = 6.0
    C2 = 7.5
    G = 2.5
    numerator = nir_band_float - red_band_float
    denominator = nir_band_float + C1 * red_band_float - C2 * blue_band_float + L
    evi = np.where(denominator == 0, 0, G * (numerator / denominator))
    return evi

def calcular_savi(nir_band, red_band):
    """
    Calcula el Índice de Vegetación Ajustado al Suelo (SAVI).
    """
    nir_band_float = nir_band.astype(np.float32)
    red_band_float = red_band.astype(np.float32)
    L = 0.5
    numerator = nir_band_float - red_band_float
    denominator = nir_band_float + red_band_float + L
    savi = np.where(denominator == 0, 0, (numerator / denominator) * (1 + L))
    return savi

def calcular_gndvi(nir_band, green_band):
    """
    Calcula el Índice de Vegetación de Diferencia Normalizada Verde (GNDVI).
    """
    nir_band_float = nir_band.astype(np.float32)
    green_band_float = green_band.astype(np.float32)
    numerator = nir_band_float - green_band_float
    denominator = nir_band_float + green_band_float
    gndvi = np.where(denominator == 0, 0, numerator / denominator)
    return gndvi

def calcular_ndre(nir_band, vre_band):
    """
    Calcula el Índice de Vegetación de Diferencia Normalizada del Borde Rojo (NDRE).
    """
    nir_band_float = nir_band.astype(np.float32)
    vre_band_float = vre_band.astype(np.float32)
    numerator = nir_band_float - vre_band_float
    denominator = nir_band_float + vre_band_float
    ndre = np.where(denominator == 0, 0, numerator / denominator)
    return ndre

def calcular_vci(ndvi_actual, ndvi_min_historico, ndvi_max_historico):
    """
    Calcula el Índice de Condición de la Vegetación (VCI).
    """
    ndvi_actual_float = ndvi_actual.astype(np.float32)
    ndvi_min_historico_float = np.asarray(ndvi_min_historico).astype(np.float32)
    ndvi_max_historico_float = np.asarray(ndvi_max_historico).astype(np.float32)
    numerator = ndvi_actual_float - ndvi_min_historico_float
    denominator = ndvi_max_historico_float - ndvi_min_historico_float
    vci = np.where(denominator == 0, 0, (numerator / denominator) * 100)
    vci = np.clip(vci, 0, 100)
    return vci

# --- Funciones Auxiliares ---
def clip_and_reproject(src_file_path, geojson_gdf, dst_crs):
    """
    Recorta una banda raster con un polígono GeoJSON y la reproyecta.
    Devuelve la banda reproyectada como un array NumPy, sus metadatos y la ruta del archivo de salida.
    """
    print(f"Procesando Banda {os.path.basename(src_file_path).split('_B')[1].split('_')[0]} desde: {src_file_path}")
    output_dir = r'C:/Users/Familia/Desktop/QGIS 2024/Salida_Indices'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with rasterio.open(src_file_path) as src:
        gdf_reprojected = geojson_gdf.to_crs(src.crs)
        out_image, out_transform = mask(src, [geom for geom in gdf_reprojected.geometry], crop=True)
        
        if out_image.shape[1] == 0 or out_image.shape[2] == 0:
            raise ValueError("El polígono no se superpone con la imagen raster.")
        
        band_dtype = getattr(src, 'dtype', None)
        band_count = getattr(src, 'count', None)

        clipped_meta = src.meta.copy()
        clipped_meta.update({
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform,
            "dtype": rasterio.float32,
            "count": band_count if band_count is not None else 1,
            "driver": "GTiff"
        })

    final_transform, final_width, final_height = calculate_default_transform(
        clipped_meta['crs'], dst_crs, clipped_meta['width'], clipped_meta['height'],
        *rasterio.transform.array_bounds(clipped_meta['height'], clipped_meta['width'], clipped_meta['transform'])
    )

    final_image = np.zeros((1, final_height, final_width), dtype=rasterio.float32)
    
    reproject(
        source=out_image,
        destination=final_image,
        src_transform=clipped_meta['transform'],
        src_crs=clipped_meta['crs'],
        dst_transform=final_transform,
        dst_crs=dst_crs,
        resampling=Resampling.bilinear
    )

    final_meta = clipped_meta.copy()
    final_meta.update({
        'crs': dst_crs,
        'transform': final_transform,
        'width': final_width,
        'height': final_height,
        'count': 1
    })
    
    final_filename = os.path.basename(src_file_path).replace('.jp2', '_reprojected.tif')
    final_path = os.path.join(output_dir, final_filename)
    
    with rasterio.open(final_path, 'w', **final_meta) as dst_final:
        dst_final.write(final_image)
        
    return final_image, final_meta, final_path

def resample_band(band_path, target_resolution, target_transform, target_crs, target_shape, dst_crs):
    """
    Remuestrea una banda raster a una resolución y extensión objetivo.
    """
    print(f"Remuestreando {os.path.basename(band_path)} a {target_resolution}m...")
    
    if not os.path.exists(band_path):
        print(f"ERROR: El archivo '{band_path}' no existe. No se puede remuestrear.")
        raise FileNotFoundError(f"Archivo no encontrado para remuestreo: {band_path}")

    with rasterio.open(band_path) as src:
        reprojected_data = np.zeros(target_shape, dtype=rasterio.float32)

        reproject(
            source=rasterio.band(src, 1),
            destination=reprojected_data,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=target_transform,
            dst_crs=dst_crs,
            resampling=Resampling.bilinear,
            num_threads=os.cpu_count()
        )
        return reprojected_data

def stretch_contrast(band_data, lower_percentile=1, upper_percentile=99):
    """
    Aplica un estiramiento de contraste lineal a una banda.
    """
    valid_data = band_data[band_data != 0]
    if valid_data.size == 0:
        return band_data, 0, 0

    min_val = np.percentile(valid_data, lower_percentile)
    max_val = np.percentile(valid_data, upper_percentile)

    stretched_band = np.interp(band_data, (min_val, max_val), (0, 255)).astype(np.uint8)
    return stretched_band, min_val, max_val

def plot_grayscale_bands(bands_dict, geojson_gdf):
    """
    Visualiza las 6 bandas en un único plot, cada una en escala de grises.
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    band_order = ['B02', 'B03', 'B04', 'B05', 'B08', 'B11']
    band_names = {
        'B02': 'Azul (B02)',
        'B03': 'Verde (B03)',
        'B04': 'Rojo (B04)',
        'B05': 'Borde Rojo (B05)',
        'B08': 'Infrarrojo Cercano (NIR) (B08)',
        'B11': 'Infrarrojo de Onda Corta (SWIR) (B11)'
    }
    
    print("\n--- Generando Plot de las 6 Bandas en Escala de Grises ---")

    for i, band_id in enumerate(band_order):
        if band_id not in bands_dict:
            print(f"Error: Banda {band_id} no encontrada. Saltando plot.")
            continue
            
        band_data = bands_dict[band_id]['data'].squeeze()
        meta = bands_dict[band_id]['meta']
        
        stretched_data, _, _ = stretch_contrast(band_data)
        data_to_plot = np.ma.masked_equal(stretched_data, 0)
        
        transform = meta['transform']
        extent = (transform.c, transform.c + stretched_data.shape[1] * transform.a,
                  transform.f + stretched_data.shape[0] * transform.e, transform.f)
        
        im = axes[i].imshow(data_to_plot, cmap='gray', extent=extent)
        
        axes[i].set_title(band_names[band_id])
        axes[i].set_xlabel('Píxeles (Este)')
        axes[i].set_ylabel('Píxeles (Norte)')
        
        axes[i].set_facecolor('dimgray')
        
        try:
            gdf_reprojected = geojson_gdf.to_crs(meta['crs'])
            gdf_reprojected.plot(ax=axes[i], facecolor='none', edgecolor='lime', linewidth=2, linestyle='-')
        except Exception as e:
            print(f"Advertencia: No se pudo plotear el polígono sobre el mapa de la banda {band_id}. Error: {e}")
        
    plt.tight_layout()
    plt.show()

def process_and_plot_index(index_data, meta, index_type, output_dir, geojson_gdf):
    """
    Guarda un índice calculado como un archivo GeoTIFF y lo visualiza con el polígono superpuesto.
    """
    print(f"\n--- Procesando y ploteando el índice: {index_type} ---")
    
    # Definición de mapas de color y rangos específicos para cada índice
    cmap = None
    vmin, vmax = None, None

    if index_type == 'NDVI':
        cmap = 'RdYlGn'
        vmin, vmax = -1, 1
    elif index_type == 'NDWI':
        cmap = 'GnBu'
        vmin, vmax = -1, 1
    elif index_type == 'NDMI':
        colors_ndmi = ["red", "darkorange", "gold", "yellowgreen", "forestgreen"]
        nodes_ndmi_original = [-1.0, -0.2, 0.0, 0.4, 1.0]
        nodes_ndmi_mapped = [(x + 1) / 2 for x in nodes_ndmi_original]
        cmap = LinearSegmentedColormap.from_list("custom_ndmi_cmap", list(zip(nodes_ndmi_mapped, colors_ndmi)))
        vmin, vmax = -1, 1
    elif index_type == 'EVI':
        cmap = 'YlGn'
        vmin, vmax = 0, 1
    elif index_type == 'SAVI':
        cmap = 'RdYlGn'
        vmin, vmax = -1, 1
    elif index_type == 'GNDVI':
        cmap = 'YlGn'
        vmin, vmax = -1, 1
    elif index_type == 'NDRE':
        cmap = 'YlGn'
        vmin, vmax = -1, 1
    elif index_type == 'VCI':
        colors_vci = ["red", "orange", "gold", "yellowgreen", "forestgreen", "darkgreen"]
        nodes_vci = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        cmap = LinearSegmentedColormap.from_list("custom_vci_cmap", list(zip(nodes_vci, colors_vci)))
        vmin, vmax = 0, 100
    else:
        print(f"Advertencia: Tipo de índice '{index_type}' no tiene mapa de color predefinido.")
        cmap = 'viridis'
        vmin, vmax = 0, 1
    
    # 1. Guardar el archivo GeoTIFF
    index_data = index_data.astype(np.float32)
    index_meta = meta.copy()
    index_meta.update({
        'dtype': 'float32',
        'count': 1,
        'nodata': 0
    })

    index_filename = os.path.join(output_dir, f'{index_type}.tif')
    with rasterio.open(index_filename, 'w', **index_meta) as dst:
        dst.write(index_data, 1)
    
    print(f"Índice {index_type} guardado en: {index_filename}")

    # 2. Generar el plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    transform = meta['transform']
    extent = (transform.c, transform.c + index_data.shape[1] * transform.a,
              transform.f + index_data.shape[0] * transform.e, transform.f)
    
    data_to_plot = np.ma.masked_equal(index_data, 0)
    im = ax.imshow(data_to_plot, cmap=cmap, vmin=vmin, vmax=vmax, extent=extent)
    
    ax.set_title(f'Mapa de {index_type} del Lote Agrícola')
    ax.set_xlabel('Píxeles (Este)')
    ax.set_ylabel('Píxeles (Norte)')
    ax.set_facecolor('black') # Fondo contrastante
    fig.colorbar(im, ax=ax, label=f'Valor {index_type}')

    try:
        gdf_reprojected = geojson_gdf.to_crs(meta['crs'])
        gdf_reprojected.plot(ax=ax, facecolor='none', edgecolor='lime', linewidth=2) # Línea fosforescente
    except Exception as e:
        print(f"Advertencia: No se pudo plotear el polígono sobre el mapa del índice {index_type}. Error: {e}")
    
    plt.tight_layout()
    plt.show()


# --- Configuración de Rutas y Nombres de Archivo ---
base_dir_all_bands = r'C:\Users\Familia\Desktop\QGIS 2024\20240613\20240613B2-3-4-5-8y11'
geojson_path = r'C:\Users\Familia\ELIO PYTHON FINAL\lote_agricola_20240613_1.geojson' # <--- ¡Modifica el nombre aquí con el de tu archivo!
output_dir = r'C:/Users/Familia/Desktop/QGIS 2024/Salida_Indices'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

dst_crs = 'EPSG:4326'

# --- Script Principal ---
if __name__ == "__main__":
    try:
        geojson_gdf = load_geojson(geojson_path)
        if geojson_gdf is None:
            raise Exception("No se pudo cargar el archivo GeoJSON.")

        band_files = {
            'B02': 'T20JKT_20240613T142719_B02_10m.jp2', # Azul
            'B03': 'T20JKT_20240613T142719_B03_10m.jp2', # Verde
            'B04': 'T20JKT_20240613T142719_B04_10m.jp2', # Rojo
            'B05': 'T20JKT_20240613T142719_B05_20m.jp2', # Red Edge (20m)
            'B08': 'T20JKT_20240613T142719_B08_10m.jp2', # NIR
            'B11': 'T20JKT_20240613T142719_B11_20m.jp2'  # SWIR1 (20m)
        }

        bandas_reproyectadas = {}
        
        print("\n--- Iniciando Reproyección de Bandas a WGS84 ---")
        for band_id, filename in band_files.items():
            ruta_banda = os.path.join(base_dir_all_bands, filename)
            banda_data, banda_meta, banda_ruta_final = clip_and_reproject(ruta_banda, geojson_gdf, dst_crs)
            
            bandas_reproyectadas[band_id] = {
                'data': banda_data,
                'meta': banda_meta,
                'path': banda_ruta_final
            }
        
        print("\n--- Recorte y Reproyección de bandas completado ---")

        print("\n--- Remuestreando bandas de 20m a 10m ---")
        target_meta = bandas_reproyectadas['B02']['meta']
        target_resolution = 10
        target_transform = target_meta['transform']
        target_crs = target_meta['crs']
        target_shape = bandas_reproyectadas['B02']['data'].shape

        for band_id in ['B05', 'B11']:
            remuestreo_data = resample_band(
                bandas_reproyectadas[band_id]['path'],
                target_resolution,
                target_transform,
                target_crs,
                (1, target_shape[1], target_shape[2]),
                dst_crs
            )
            if remuestreo_data is not None:
                bandas_reproyectadas[band_id]['data'] = remuestreo_data
                bandas_reproyectadas[band_id]['meta'].update({
                    'width': target_shape[2],
                    'height': target_shape[1],
                    'transform': target_transform
                })
        print("--- Remuestreo de bandas completado ---")

        # Extraer las bandas finales (2D arrays)
        band_blue = bandas_reproyectadas['B02']['data'][0]
        band_green = bandas_reproyectadas['B03']['data'][0]
        band_red = bandas_reproyectadas['B04']['data'][0]
        band_vre = bandas_reproyectadas['B05']['data'][0]
        band_nir = bandas_reproyectadas['B08']['data'][0]
        band_swir1 = bandas_reproyectadas['B11']['data'][0]

        # --- Visualización de Bandas Individuales en Escala de Grises ---
        bands_to_plot = {
            'B02': {'data': band_blue, 'meta': bandas_reproyectadas['B02']['meta']},
            'B03': {'data': band_green, 'meta': bandas_reproyectadas['B03']['meta']},
            'B04': {'data': band_red, 'meta': bandas_reproyectadas['B04']['meta']},
            'B05': {'data': band_vre, 'meta': bandas_reproyectadas['B05']['meta']},
            'B08': {'data': band_nir, 'meta': bandas_reproyectadas['B08']['meta']},
            'B11': {'data': band_swir1, 'meta': bandas_reproyectadas['B11']['meta']}
        }
        plot_grayscale_bands(bands_to_plot, geojson_gdf)

        # --- Cálculo y Visualización de Índices ---
        meta_10m = bandas_reproyectadas['B02']['meta']
        
        print("\n--- Calculando y Visualizando Índices ---")

        # NDWI (verde y NIR)
        ndwi_result = calcular_ndwi(band_green, band_nir)
        process_and_plot_index(ndwi_result, meta_10m, 'NDWI', output_dir, geojson_gdf)
        
        # NDVI (NIR y rojo)
        ndvi_result = calcular_ndvi(band_nir, band_red)
        process_and_plot_index(ndvi_result, meta_10m, 'NDVI', output_dir, geojson_gdf)
        
        # NDMI (NIR y SWIR1)
        ndmi_result = calcular_ndmi(band_nir, band_swir1)
        process_and_plot_index(ndmi_result, meta_10m, 'NDMI', output_dir, geojson_gdf)
        
        # EVI (NIR, rojo y azul)
        evi_result = calcular_evi(band_nir, band_red, band_blue)
        process_and_plot_index(evi_result, meta_10m, 'EVI', output_dir, geojson_gdf)
        
        # SAVI (NIR y rojo)
        savi_result = calcular_savi(band_nir, band_red)
        process_and_plot_index(savi_result, meta_10m, 'SAVI', output_dir, geojson_gdf)
        
        # GNDVI (NIR y verde)
        gndvi_result = calcular_gndvi(band_nir, band_green)
        process_and_plot_index(gndvi_result, meta_10m, 'GNDVI', output_dir, geojson_gdf)
        
        # NDRE (NIR y borde rojo)
        ndre_result = calcular_ndre(band_nir, band_vre)
        process_and_plot_index(ndre_result, meta_10m, 'NDRE', output_dir, geojson_gdf)

        # VCI (usa el NDVI)
        ndvi_values_for_percentiles = ndvi_result[ndvi_result != 0]
        if ndvi_values_for_percentiles.size > 0:
            ndvi_min_simulado = np.percentile(ndvi_values_for_percentiles, 2)
            ndvi_max_simulado = np.percentile(ndvi_values_for_percentiles, 98)
        else:
            ndvi_min_simulado = 0.0
            ndvi_max_simulado = 1.0
            print("Advertencia: No hay píxeles de NDVI válidos para VCI. Usando 0-1 como min/max simulado.")
        
        vci_result = calcular_vci(ndvi_result, ndvi_min_simulado, ndvi_max_simulado)
        process_and_plot_index(vci_result, meta_10m, 'VCI', output_dir, geojson_gdf)
        
        print("--- Cálculo y Visualización de Índices completado ---")

    except Exception as e:
        print(f"\nOcurrió un error inesperado: {e}")
        print("Asegúrate de que los archivos GeoJSON y de las bandas no estén dañados.")
    finally:
        print("\n--- Iniciando limpieza de archivos .tif generados ---")
        suffixes_to_delete = ('_reprojected.tif', '.tif')
        for item in os.listdir(output_dir):
            if item.endswith(suffixes_to_delete):
                file_to_delete = os.path.join(output_dir, item)
                try:
                    if os.path.exists(file_to_delete):
                        os.remove(file_to_delete)
                        print(f"Archivo .tif generado eliminado: {file_to_delete}")
                except Exception as e:
                    print(f"Error al eliminar el archivo {file_to_delete}: {e}")
        print("--- Limpieza de archivos .tif generados completada ---")