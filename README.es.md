[English](README.md) | **Español**

# Modelos de Priorización para la Conservación

## Resumen

Este repositorio contiene scripts para preparar datos espaciales y ejecutar análisis sistemáticos de priorización para la conservación en cuatro regiones geográficas de Colombia: dos análisis de escala nacional (terrestre y marino) y dos análisis regionales en paisajes prioritarios de conservación (SIRAP Eje Cafetero y SIRAP Orinoquía) — cubriendo el flujo de trabajo completo, desde la preparación de datos crudos hasta los resultados de priorización ya resueltos. Todos los modelos se construyen con el paquete de R [`prioritizr`](https://prioritizr.net/), formulando la planificación espacial de la conservación como programas lineales enteros mixtos (MILP, por sus siglas en inglés). Los problemas se resuelven con Gurobi por defecto, aunque también se admiten solvers de código abierto (ver [Requisitos](#requisitos)).

**Regiones evaluadas actualmente:**

| Región | Ámbito | Escala | Resolución | Script |
|----|----|----|----|----|
| Colombia (nacional) | Terrestre | Nacional | 1 km | `2_1_national_terrestrial_model.R` |
| Colombia (nacional) | Marino | Nacional | 1 km | `2_2_national_marine_model.R` |
| SIRAP Eje Cafetero | Terrestre | Regional | 300 m | `3_1_sirap_eje_cafetero_model.R` |
| SIRAP Orinoquía | Terrestre | Regional | 500 m | `3_2_sirap_orinoquia_model.R` |

Se construyó un modelo independiente para cada una de las cuatro regiones, utilizando distintos datos de entrada y generando soluciones de diferente resolución. El objetivo es crear modelos para cada uno de los SIRAP territoriales y temáticos, utilizando datos y parámetros de modelación relevantes a nivel local; actualmente solo se han creado dos.

<!-- TODO: agregar mapas de cada región aquí para visualizar las áreas -->

## Estructura del repositorio

```
.
└── scripts/
    ├── utils.R                              # Funciones auxiliares compartidas, plantillas de unidades de planificación y tablas de escenarios que usan todos los scripts
    ├── 1_data_prep.R                        # Limpieza, recorte y preparación de capas (ejecutar primero, requerido por todos los modelos)
    ├── 2_1_national_terrestrial_model.R     # Modelo de priorización nacional terrestre
    ├── 2_2_national_marine_model.R          # Modelo de priorización nacional marino
    ├── 3_1_sirap_eje_cafetero_model.R       # Modelo regional: SIRAP Eje Cafetero
    └── 3_2_sirap_orinoquia_model.R          # Modelo regional: SIRAP Orinoquía
```

Todos los scripts asumen que el directorio de trabajo es la **raíz del repositorio** (usan `source("scripts/utils.R")`), y `utils.R` a su vez construye varias otras rutas a partir de esa raíz. Ejecuta los scripts desde un proyecto de R/`.Rproj` ubicado en la raíz, o usa `setwd()` para posicionarte ahí antes de correrlos.

<!-- TODO: si terminan agregando 1 o 2 scripts más, incluirlos aquí (p. ej. un script de post-procesamiento o visualización) -->

## Datos

Los datos de entrada **no están incluidos en este repositorio** debido a su tamaño. Todos los conjuntos de datos deben descargarse o solicitarse por separado y colocarse en la estructura de directorios local esperada (más abajo) antes de ejecutar `1_data_prep.R`.

| Conjunto de datos | Uso | Fuente / notas |
|----|----|----|
| IHEH (Índice de Huella Humana) 2022 y proyección 2030 | Capa de costo (terrestre y SIRAP) | [Producido por el Instituto Humboldt](https://geonetwork.humboldt.org.co/geonetwork/srv/spa/catalog.search#/metadata/7d8f0aeb-8136-45a7-a469-f0016f618250) (IAVH, 2024. Contrato de prestación de servicios No. 23-22/187-23/0017-197PS. Instituto Humboldt. Bogotá, D. C., Colombia.) |
| IHEH 2030 | Capa de costo (terrestre y SIRAP) | [Producido por el Instituto Humboldt](https://reporte.humboldt.org.co/biodiversidad/2019/cap2/203/#seccion10); datos proporcionados directamente por los autores del informe. |
| Huella humana marina | Capa de costo (marino) | Elaborada por INVEMAR y revisada por la Mesa Nacional. Esta capa combina presiones provenientes de factores que contribuyen a la pérdida de biodiversidad en las zonas marinas y costeras del Caribe y el Pacífico colombianos, usando datos de línea base 2019–2020 a 300 m de resolución. Considera nueve variables: 1) transformación de hábitats costeros, 2) acceso humano a la zona costera, 3) ocupación humana, 4) luz nocturna, 5) contaminación marina, 6) pesca industrial, 7) pesca artesanal, 8) tráfico marítimo, y 9) actividad portuaria. |
| RUNAP (Registro Único Nacional de Áreas Protegidas) | Restricción de áreas ya incluidas (locked-in, todas las regiones) | [Registro Único Nacional AP](https://www.datos.gov.co/dataset/runap-Registro-Unico-Nacional-AP/n9kx-xwgg/about_data) (última actualización 9 de abril de 2026) |
| WDOECM — Otras Medidas Efectivas de Conservación (OMEC) | Restricción de áreas ya incluidas (todas las regiones) | [Protected Planet / WDOECM](https://www.protectedplanet.net/en/thematic-areas/oecms), filtrado a Colombia (descargado en junio de 2026) |
| Ecosistemas marinos | Capa de features (marino) | Consolidada por INVEMAR, combinando sus análisis de unidades de ecosistemas de zonas profundas y someras en una sola capa; no alojada públicamente. |
| Manglares | Capa de features (marino) | [Producido por INVEMAR](https://acceso-datos-ambientales-invemar.hub.arcgis.com/datasets/INVEMAR::manglares-colombia/about) (última actualización 22 mar. 2023) |
| Ecosistemas continentales y costeros (biomas IAvH) 2024 | Capa de features (terrestre) | [Producido por IDEAM](https://www.ideam.gov.co/ecosistemas) — *Mapa de Ecosistemas Continentales, Costeros y Marinos (MEC)*, 100K, 2024 |
| Páramos | Capa de features (terrestre y SIRAP) | Producido por el MADS y [obtenido de SIAC](https://siac-datosabiertos-mads.hub.arcgis.com/maps/a85c9d818c9640649318a33b9fa115fd/about) (última actualización 27 ago. 2026) |
| Bosque seco tropical | Capa de features (terrestre y SIRAP) | Producido por el MADS y [obtenido de SIAC](https://siac-datosabiertos-mads.hub.arcgis.com/maps/a85c9d818c9640649318a33b9fa115fd/about) (última actualización 27 ago. 2026) |
| Humedales (nacional) | Capa de features (terrestre y SIRAP) | Producido por el MADS y [obtenido de SIAC](https://siac-datosabiertos-mads.hub.arcgis.com/maps/a85c9d818c9640649318a33b9fa115fd/about) (última actualización 27 ago. 2026) |
| Humedales (específico de Eje Cafetero) | Capa de features (Eje Cafetero) | Datos complementarios de humedales dentro de Eje Cafetero, proporcionados directamente, no alojados públicamente |
| Sabanas (SIRAP Orinoquía) | Capa de features (Orinoquía) | Datos complementarios de sabanas en Orinoquía, proporcionados directamente, no alojados públicamente |
| Congriales (SIRAP Orinoquía) | Capa de features (Orinoquía) | Datos proporcionados directamente, no alojados públicamente |
| Carbono (AGB + BGB) | Capa de features (terrestre) | Spawn et al. 2020 — [conjunto de datos](https://doi.org/10.3334/ORNLDAAC/1763); el repositorio usa actualmente una versión preprocesada por un colaborador ("Jaime") — **falta documentar el método de ese preprocesamiento** |
| Zonas de recarga de agua subterránea | Capa de features (terrestre) | [Fuente: IDEAM](https://bart.ideam.gov.co/cneideam/Capasgeo/#:~:text=Zonas%5FPotenciales%5Fde%5FRecarga%5Fde%5FAgua%5FSubterraneas%5FEna2018), *Zonas Potenciales de Recarga de Aguas Subterráneas*, ENA 2018 |
| Cobertura de la tierra (IDEAM CLC 2022) | Referencia/visualización | IDEAM, Cobertura de la Tierra 100K, 2022 — [portal de datos abiertos ArcGIS](https://experience.arcgis.com/experience/568ddab184334f6b81a04d2fe9aac262/page/Datos-Abiertos-Geogr%C3%A1ficos-/) |
| Modelos de distribución de especies (BioModelos) | Capa de features (terrestre) | [Conjunto de datos BioModelos](https://geonetwork.humboldt.org.co/geonetwork/srv/spa/catalog.search#/metadata/0a1a6bdf-3231-4a77-8031-0dc3fa40f21b), alojado por el IAvH |
| Evaluaciones de la Lista Roja de la UICN | Capa de features de especies (terrestre) y definición de objetivos de conservación | Descarga masiva de la [Lista Roja de la UICN](https://www.iucnredlist.org/). Se usó el estado de amenaza determinado a nivel nacional cuando estaba disponible, y los vacíos se llenaron con la Lista Roja de la UICN. |
| Tierras de comunidades afrocolombianas | Referencia/visualización | [Producido por la ANT](https://data-agenciadetierras.opendata.arcgis.com/datasets/agenciadetierras::consejo-comunitario-titulado/about) (Agencia Nacional de Tierras) |
| Resguardos indígenas | Referencia/visualización | [Producido por la ANT](https://data-agenciadetierras.opendata.arcgis.com/datasets/agenciadetierras::resguardo-indigena-formalizado/about) (Agencia Nacional de Tierras) |
| Sitios RAMSAR | Referencia/visualización | [Sitios RAMSAR de Colombia](https://siac-datosabiertos-mads.hub.arcgis.com/maps/3b63a187479543eb858028ecaa9d068b/about), producido por el MADS |
| Reservas de Biosfera UNESCO | Referencia/visualización | [Reservas de Biosfera de Colombia](https://siac-datosabiertos-mads.hub.arcgis.com/maps/3b63a187479543eb858028ecaa9d068b/about), producido por el MADS |
| Reservas Forestales (Ley 2ª de 1959) | Referencia/visualización | [Datos producidos por el MADS](https://siac-datosabiertos-mads.hub.arcgis.com/maps/3b63a187479543eb858028ecaa9d068b/about) |
| Zonas de Reserva Campesina Constituida | Referencia/visualización | [Producido por la ANT](https://data-agenciadetierras.opendata.arcgis.com/datasets/zonas-de-reserva-campesina-constituida/explore?location=3.091663%2C-71.907458%2C6) |
| ECC Eje Cafetero (*Estructura Ecológica Complementaria*) | Referencia/visualización (Eje Cafetero) | Datos proporcionados directamente, no alojados públicamente |
| Definiciones de escenarios y objetivos del modelo (`corridas_*.xlsx`) | Define los features de conservación, objetivos, capas de áreas ya incluidas y el costo usado en cada corrida del modelo | Hojas de cálculo proporcionadas por la Mesa Nacional — una por región: `corridas_05062026.xlsx` (nacional terrestre), `corridas_SIRAP_EC_16072026.xlsx` (Eje Cafetero), `corridas_SIRAP_ORI_16072026.xlsx` (Orinoquía). Los escenarios nacionales marinos se construyen manualmente en `utils.R`. Son **insumos requeridos**, no datos SIG crudos — sin ellos los scripts de modelo no tienen sobre qué iterar. <!-- TODO: pendiente el permiso de los socios, estos archivos se subirán directamente al repositorio (p. ej. bajo data/model_inputs/) para que vengan incluidos con `git clone` — actualizar esta fila y la sección de Configuración inicial cuando eso suceda --> |
| AICAs / KBAs | Feature planeada, aún no implementada | Acceso a los datos pendiente en la versión actual de los scripts |

### Estructura de directorio local esperada {#estructura-de-directorio-local-esperada}

Aunque no forma parte del repositorio, los scripts asumen esta organización relativa a la raíz del proyecto (las rutas provienen directamente de `utils.R` / `1_data_prep.R`):

```
data/
├── costs/                        # Rásteres IHEH (2022 y 2030), huella humana marina/
├── includes/                     # Shapefiles de RUNAP, OMEC, comunidades y resguardos indígenas
├── features/                     # Ecosistemas, páramos, bosque seco, humedales, manglares, carbono, recarga de agua, capas específicas de Orinoquía
├── sirap_actualizado/            # Shapefiles de límites de SIRAP/territoriales (usados para construir las plantillas de unidades de planificación de SIRAP)
├── visualization/                # Capas usadas solo para visualización: ECC, ZRC, RAMSAR, reservas de biosfera, reservas forestales, IDEAM CLC
├── redlist_species_data_.../     # Descarga masiva de la Lista Roja de la UICN (assessments.csv)
├── temp_outputs/                 # Resultados intermedios escritos por los scripts (se crea automáticamente)
│   ├── national/
│   └── sirap/{eje_cafetero,orinoquia}/
├── model_input_lyrs/             # GeoTIFFs/shapefiles preparados por región (se crea automáticamente). Se usa para visualizar los insumos preparados del modelo.
│   ├── national/
│   └── sirap/{eje_cafetero,orinoquia}/
└── model_inputs/                 # Matrices (.rds) + hojas de cálculo de escenarios usadas directamente en los modelos de prioritizr (se crea automáticamente, excepto los archivos Excel)
    ├── national/                 #   -> también contiene corridas_05062026.xlsx
    └── sirap/
        ├── eje_cafetero/         #   -> también contiene corridas_SIRAP_EC_16072026.xlsx
        └── orinoquia/            #   -> también contiene corridas_SIRAP_ORI_16072026.xlsx

results/                          # Resultados de los modelos (se crea automáticamente en cada script de modelo)
├── national/{terrestrial,marine}/
└── sirap/{eje_cafetero,orinoquia}/
```

`temp_outputs/`, `model_input_lyrs/`, `model_inputs/` y `results/...` se crean automáticamente si no existen — solo necesitas poblar de antemano `costs/`, `includes/`, `features/`, `sirap_actualizado/`, `visualization/`, `redlist_species_data_.../`, y los tres archivos de escenarios `corridas_*.xlsx`.

## Requisitos {#requisitos}

- **R** con el paquete [`pacman`](https://cran.r-project.org/package=pacman) (se instala automáticamente si falta; se usa para instalar el resto de paquetes automáticamente).
- **Paquetes clave de R** (instalados automáticamente mediante `pacman::p_load()` donde es posible): `tidyverse`, `readxl`, `terra`, `sf`, `purrr`, `Matrix`, `prioritizr`, `exactextractr`.
- **Gurobi** (con licencia válida) — los scripts de modelo actualmente llaman a `add_gurobi_solver()`. Es un solver comercial: necesitarás una licencia (Gurobi ofrece [licencias académicas gratuitas](https://www.gurobi.com/academia/academic-program-and-licenses/)) y el paquete de R `gurobi` instalado desde la instalación de Gurobi misma (no está en CRAN, así que `pacman`/`p_load` no puede instalarlo automáticamente — instálalo manualmente siguiendo las instrucciones de configuración de R de Gurobi). `prioritizr` también admite solvers de código abierto (p. ej. `Rsymphony`, `HiGHS` mediante `highs`, o `cbc` mediante `rcbc`) — puedes reemplazar las llamadas a `add_gurobi_solver()` en los cuatro scripts de modelo por alguno de estos si prefieres evitar el requisito de licencia.

## Configuración inicial

1. Clona este repositorio y ábrelo como directorio de trabajo de R (p. ej. mediante un `.Rproj` en la raíz).
2. Crea las subcarpetas de `data/` indicadas en [Estructura de directorio local esperada](#estructura-de-directorio-local-esperada) y llénalas con los conjuntos de datos de [Datos](#datos), incluyendo los tres archivos de escenarios `corridas_*.xlsx`.
3. Instala Gurobi y su paquete de R manualmente (si lo vas a usar; ver [Requisitos](#requisitos)) — es la única dependencia que los scripts no pueden instalar por ti.
4. Ejecuta `scripts/1_data_prep.R`. Esto es obligatorio antes de ejecutar **cualquier** script de modelo — construye las plantillas de unidades de planificación y todas las matrices de costo/features preparadas que usan las cuatro regiones. Nota: este script, en particular la sección de rangos de especies de BioModelos, puede tardar **muchas horas (8 o más)**; considera ejecutarlo durante la noche o en otra máquina.

## Uso / orden de ejecución

```
1_data_prep.R                       # Obligatorio primero — prepara los datos de todas las regiones
   │
   ├── 2_1_national_terrestrial_model.R
   ├── 2_2_national_marine_model.R
   ├── 3_1_sirap_eje_cafetero_model.R
   └── 3_2_sirap_orinoquia_model.R
```

Una vez ejecutado `1_data_prep.R`, los cuatro scripts de modelo son independientes entre sí y pueden ejecutarse en cualquier orden, individualmente o todos juntos, según la(s) región(es) que necesites. `utils.R` se carga automáticamente desde los demás scripts y no necesita ejecutarse por separado.

Cada script de modelo:

1. Carga `utils.R` y los insumos preparados de la región, junto con su tabla de escenarios (proveniente del `corridas_*.xlsx` correspondiente, ya convertida en data frame dentro de `utils.R`).
2. Itera sobre cada escenario (fila) de esa tabla con `purrr::pmap()`, construyendo y resolviendo un problema de `prioritizr` para cada uno — fijando las capas indicadas como ya incluidas (p. ej. RUNAP, RUNAP+OMEC), aplicando los objetivos de features especificados, y usando la capa de costo especificada.
3. Resuelve con Gurobi (`gap = 0.05`) y escribe un ráster con la solución y estadísticas resumen para cada escenario.
4. **Las corridas se pueden reanudar**: al inicio de la sección de ejecución de cada script, se revisan `master_eval_summary.csv` y `failed_scenarios.txt` de la región, y se omite cualquier escenario ya completado o ya registrado como fallo (que no sea de presolve). Un segundo paso reintenta los escenarios que fallaron solo en la verificación de presolve, forzando la resolución.

## Resultados

Se escriben en `results/<región>/` (se crea automáticamente), una carpeta por región:

| Archivo | Contenido |
|----|----|
| `<model_name>.tif` | Solución rasterizada de ese escenario: celdas nuevas seleccionadas como prioritarias vs. celdas ya incluidas (áreas ya protegidas). Las celdas que no forman parte de la solución quedan como NA. |
| `<model_name>_summary.csv` | Detalle de cobertura de objetivos por escenario, listando cada feature como una fila (incluyendo especies y ecosistemas individuales, según aplique) con detalles de cuánto queda representado en la solución. |
| `master_eval_summary.csv` | Una fila por escenario en toda la región: conteo de celdas totales/nuevas/ya incluidas, costo, % de objetivos cumplidos. Se agrega (no se sobrescribe) conforme los escenarios se completan, lo cual permite reanudar las corridas. |
| `failed_scenarios.txt` | Registro de nombre de escenario + marca de tiempo + mensaje de error para cualquier escenario que no se pudo resolver. Se usa para reanudar el modelo si se interrumpe a medio camino. |
| `resultados_todos.csv` *(solo scripts de SIRAP)* | Todos los CSV de resumen por escenario de la región, concatenados en un solo archivo |

`model_name` codifica los objetivos y configuración del escenario (p. ej. objetivos de features, qué capas de áreas ya incluidas se usaron, y la capa de costo empleada) — ver las funciones `build_model_name()` en `utils.R` para la lógica exacta de nomenclatura de cada región.

## Notas / advertencias

- El tiempo de ejecución de Gurobi escala con el número de escenarios y de unidades de planificación; las corridas a escala nacional tomarán considerablemente más tiempo que las regiones SIRAP. Algunos escenarios nacionales terrestres evalúan los ~8,000 objetivos de especies y pueden requerir **más de 64 GB de RAM** — se recomienda ejecutar los scripts nacionales en una máquina virtual con mayor memoria/capacidad de procesamiento, en lugar de un portátil típico. `failed_scenarios.txt` documentará qué escenarios fallaron por limitaciones de memoria.
- Varias capas de datos (sabanas y congriales de Orinoquía, humedales de Eje Cafetero, ecosistemas marinos y huella humana marina) fueron compartidas directamente por los socios y no provienen de un portal público — se indica en [Datos](#datos) donde ese sea el caso.
- La sección de especies de BioModelos actualmente solo se ejecuta para el modelo **nacional terrestre**; los modelos SIRAP aún no incluyen especies como feature, pero una evaluación posterior en esos scripts sí evalúa la cobertura a nivel de especie en las soluciones y la agrega a su `<model_name>_summary.csv`.

<!-- TODO: agregar cualquier otra cosa que valga la pena señalar — p. ej. vacíos de datos conocidos, supuestos en la definición de objetivos, a quién contactar para acceder a los conjuntos de datos compartidos/internos -->

## Cita / contacto

<!-- TODO: cómo te gustaría que se citara este trabajo, dando crédito a las personas involucradas. -->