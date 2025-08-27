# indices-espectrales-python
indices-espectrales-python: Automatiza el cálculo de índices (NDVI, EVI, SAVI, NDMI, NDWI, GNDVI, NDRE y VCI) para monitoreo agroambiental con imágenes Sentinel-2. Incluye herramientas para selección de lotes (poligonos) con GeoJSON, procesamiento de bandas y generación de mapas de los difernetes indices de vegetación y humedad. 🛰️🌱

## Cómo usar este repositorio?
### Fases del Flujo de Trabajo
El proceso para el monitoreo agroambiental con este script se divide en dos fases principales, que se corresponden con las dos partes del código.
### Fase 1: Delimitación del Área de Interés (AOI)
Esta fase es crucial para definir el lote agrícola o el área específica que se desea analizar. El script genera un mapa interactivo en un archivo HTML que permite al usuario dibujar un polígono sobre la zona de interés.

1. Configuración de la fecha: Antes de ejecutar, debes abrir el script y ajustar la variable FECHA_DE_IMAGEN. Asegúrate de que esta fecha, en formato YYYYMMDD, coincida con la de la imagen Sentinel-2 que vas a procesar. Por ejemplo, si tu imagen es del 13 de junio de 2024, la variable debe ser "20240613".

2. Generación del mapa interactivo: Al ejecutar el script, se creará un archivo HTML (mapa_lote_agricola_FECHA_DE_IMAGEN_N.html) en el mismo directorio. El sufijo _N es un número secuencial que evita que sobrescribas archivos anteriores.

3. Delineación del polígono: Abre el archivo HTML en tu navegador web. En el mapa, utiliza la herramienta de dibujo de polígonos para trazar el contorno exacto de tu lote agrícola.

4. Guardado del archivo de coordenadas: Una vez que hayas terminado de dibujar, haz clic en el botón de guardar dentro del mapa. Esto generará un archivo de coordenadas en formato GeoJSON (lote_agricola_FECHA_DE_IMAGEN_N.geojson). Este archivo es el insumo principal para la siguiente fase.

### Fase 2: Procesamiento y Cálculo de Índices Espectrales
Esta fase se activa después de que has definido tu área de interés (AOI) en la Fase 1. Se encarga de procesar las imágenes satelitales y generar los mapas de índices.

1. Preparación de los datos:
Descarga las bandas de la imagen Sentinel-2 para la fecha que seleccionaste. Las bandas necesarias son: B2, B3, B4, B5, B8, y B11.
Coloca estos archivos .jp2 en el directorio de entrada que se especifica en el script.

3. Configuración de rutas:
En el script principal, localiza las variables base_dir_all_bands y geojson_path.
Actualiza base_dir_all_bands con la ruta a la carpeta donde guardaste las bandas de Sentinel-2.
Actualiza geojson_path con la ruta completa del archivo .geojson que generaste en la Fase 1.

3. Ejecución del script:
Ejecuta el script desde tu terminal. El proceso se encargará automáticamente de reproyectar, recortar y remuestrear las bandas a una resolución de 10 metros, garantizando que todos los datos estén alineados.
Luego, calculará cada uno de los índices espectrales definidos en el código.

4. Resultados y visualización:
Durante la ejecución, el script generará varias ventanas de visualización que mostrarán los mapas de los índices calculados, con el polígono de tu lote superpuesto.
Adicionalmente, guardará una copia de cada índice en formato GeoTIFF (.tif) dentro del directorio de salida (Salida_Indices), listos para ser utilizados en software SIG como QGIS.
El script incluye una función de limpieza final que elimina los archivos .tif temporales para mantener el directorio ordenado, conservando únicamente los resultados finales.

## Autor
Elio Leguina Huertas¹, Rubén Ledesma²'³

## Afiliaciones
¹ Facultad de Ciencias Naturales - Universidad Nacional de Salta

² Facultad de Ciencias Exactas - Universidad Nacional de Salta

³ Instituto de Energía no Convencional - CONICET

## Contacto
Elio Leguina Huertas

Correo electrónico: elioleguinahuertas@gmail.com

Rubén Ledesma

Correo electrónico: rdledesma@exa.unsa.edu.ar

Este trabajo será publicado en ReTec proximamente.....
