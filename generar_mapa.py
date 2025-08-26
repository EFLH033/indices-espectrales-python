# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 01:08:23 2025

@author: Familia
"""

import folium
from folium.plugins import Draw
from datetime import datetime
import os
import glob # <--- ¡IMPORTANTE! Se agregó esta línea para importar el módulo glob

# Coordenadas aproximadas del Valle de Lerma, Salta
salta_coords = [-24.782, -65.411]

# --- ¡IMPORTANTE! Define aquí la fecha de la imagen que estás analizando ---
# Esta fecha se usará en el nombre del archivo GeoJSON y HTML.
# Asegúrate de que coincida con la fecha de tus imágenes Sentinel-2.
FECHA_DE_IMAGEN = "20240613" # Formato YYYYMMDD (Ej: 13 de junio de 2024)

# Crea el mapa base con la capa de satélite
m = folium.Map(
    location=salta_coords,
    zoom_start=11,
    tiles='Esri.WorldImagery',
    attr='Tiles &copy; Esri',
    name='Satelite'
)

# Añade una capa de etiquetas de calles y ciudades que se superpone a la de satélite
folium.TileLayer(
    'Stamen Terrain', # Esta capa tiene calles y etiquetas claras
    attr='Map tiles by Stamen Design, under CC BY 3.0. Data by OpenStreetMap, under ODbL.',
    name='Etiquetas (Calles y Ciudades)',
    control=True, # Habilita el control de capa
    overlay=True, # Define esta capa como una superposición, no como una capa base
    opacity=0.7 # Ajusta la transparencia para ver el satélite debajo
).add_to(m)

# --- Genera un nombre de archivo único con la FECHA_DE_IMAGEN y un número secuencial ---
# Busca el número secuencial más alto para los archivos GeoJSON existentes con esa fecha.
base_filename = f'lote_agricola_{FECHA_DE_IMAGEN}'
existing_files = glob.glob(f'{base_filename}_*.geojson') # Busca archivos como 'lote_agricola_20240613_*.geojson'

max_sequence = 0
for f in existing_files:
    try:
        # Extrae el número secuencial del nombre del archivo
        seq_str = f.split('_')[-1].replace('.geojson', '')
        seq_num = int(seq_str)
        if seq_num > max_sequence:
            max_sequence = seq_num
    except ValueError:
        continue # Ignora archivos que no sigan el patrón esperado

new_sequence = max_sequence + 1
geojson_filename = f'{base_filename}_{new_sequence}.geojson' # Ej: lote_agricola_20240613_1.geojson

# Agrega la herramienta de dibujo al mapa
draw = Draw(
    export=True,
    filename=geojson_filename, # <--- Usa el nombre de archivo dinámico aquí
    draw_options={
        'polyline': False,
        'marker': False,
        'circlemarker': False,
        'circle': False,
        'rectangle': False, # Desactiva el rectángulo si solo quieres polígonos irregulares
        'polygon': True # Asegura que la herramienta de polígono esté activa
    },
    edit_options={
        'edit': False,   # No permite edición de polígonos después de dibujados
        'remove': False  # No permite eliminar polígonos
    }
)
draw.add_to(m)

# Guarda el mapa interactivo en un archivo HTML
html_output_path = f'mapa_lote_agricola_{FECHA_DE_IMAGEN}_{new_sequence}.html' # Nombre de archivo HTML dinámico
m.save(html_output_path)

print(f"Se ha creado '{html_output_path}'.")
print(f"Ábrelo en tu navegador. Dibuja un polígono y al guardarlo, se creará '{geojson_filename}'.")
print(f"Para dibujar otro polígono para la fecha {FECHA_DE_IMAGEN}, ejecuta el script de nuevo para generar el siguiente archivo secuencial.")
