# indices-espectrales-python
indices-espectrales-python: Automatiza el cálculo de índices (NDVI, EVI, SAVI, NDMI, NDWI, GNDVI, NDRE y VCI) para monitoreo agroambiental con imágenes Sentinel-2. Incluye herramientas para selección de lotes (poligonos) con GeoJSON, procesamiento de bandas y generación de mapas de los difernetes indices de vegetación y humedad. 🛰️🌱

## Cómo usar este repositorio?
###Fases del Flujo de Trabajo
El proceso para el monitoreo agroambiental con este script se divide en dos fases principales, que se corresponden con las dos partes del código.
Fase 1: Delimitación del Área de Interés (AOI)
Esta fase es crucial para definir el lote agrícola o el área específica que se desea analizar. El script genera un mapa interactivo en un archivo HTML que permite al usuario dibujar un polígono sobre la zona de interés.

1.Configuración de la fecha: Antes de ejecutar, debes abrir el script y ajustar la variable FECHA_DE_IMAGEN. Asegúrate de que esta fecha, en formato YYYYMMDD, coincida con la de la imagen Sentinel-2 que vas a procesar. Por ejemplo, si tu imagen es del 13 de junio de 2024, la variable debe ser "20240613".

2.Generación del mapa interactivo: Al ejecutar el script, se creará un archivo HTML (mapa_lote_agricola_FECHA_DE_IMAGEN_N.html) en el mismo directorio. El sufijo _N es un número secuencial que evita que sobrescribas archivos anteriores.

3.Delineación del polígono: Abre el archivo HTML en tu navegador web. En el mapa, utiliza la herramienta de dibujo de polígonos para trazar el contorno exacto de tu lote agrícola.

4.Guardado del archivo de coordenadas: Una vez que hayas terminado de dibujar, haz clic en el botón de guardar dentro del mapa. Esto generará un archivo de coordenadas en formato GeoJSON (lote_agricola_FECHA_DE_IMAGEN_N.geojson). Este archivo es el insumo principal para la siguiente fase.
2. Debe enviarme su dinero

## Autor
Elio Leguina Huertas

## Contacto

Este trabajo fue publicado en ReTec ... ..... ....
