# Proyecto de Conservación de Biodiversidad en Menorca

## 📋 Descripción del Proyecto

Este proyecto implementa un modelo de optimización para la conservación de biodiversidad en la isla de Menorca. El objetivo es determinar la asignación óptima de recursos (presupuesto limitado) para expandir los hábitats de cuatro especies clave mediante:

- **Adaptación de celdas**: Modificar el hábitat de celdas para hacerlas adecuadas para especies específicas
- **Construcción de corredores**: Conectar poblaciones existentes con nuevas áreas mediante corredores ecológicos

### Especies Objetivo

1. **Atelerix algirus** (Erizo argelino) - Peso: 1.0
2. **Martes martes** (Marta) - Peso: 1.0
3. **Eliomys quercinus** (Lirón careto) - Peso: 2.0
4. **Oryctolagus cuniculus** (Conejo europeo) - Peso: 1.5

## 🏗️ Estructura del Proyecto

```text
biodiversity-conservation-menorca/
├── 1_documentation/          # Documentación del modelo y justificaciones
│   ├── model-design/         # Versiones del diseño del modelo
│   └── minimum-area-targets/ # Justificación de objetivos de área mínima
├── 2_data/                   # Datos del proyecto
│   ├── raw/                  # Datos originales (dataset.geojson)
│   └── processed/            # Datos procesados (final_dataset.csv/geojson)
├── 3_notebooks/              # Análisis exploratorio y visualizaciones
│   ├── 00_project_statement.ipynb.ipynb
│   ├── 01_eda_analysis.ipynb
│   ├── 02_visualization_maps.ipynb
│   └── iterations/           # Iteraciones del modelo (v1, v2)
├── 4_source_code/            # Código fuente principal
│   ├── data_preparation.py   # Preparación y procesamiento de datos
│   ├── optimization_model.py # Modelo de optimización CP-SAT
│   └── main.py               # Script principal de ejecución
├── 5_results/                # Resultados generados
│   ├── EDA/                  # Análisis exploratorio de datos
│   ├── solutions/            # Soluciones óptimas y mapas
│   ├── iterations/           # Resultados de iteraciones
│   └── tables/               # Tablas de resultados
└── 6_final_submission/       # Entrega final
```

## 🚀 Inicio Rápido

### Requisitos Previos

```bash
pip install pandas geopandas numpy networkx ortools matplotlib seaborn
```

O instala desde `requirements.txt` (si está disponible).

### Ejecución del Modelo

1. **Preparar los datos** (si aún no están procesados):

```bash
cd 4_source_code
python data_preparation.py
```

2. **Ejecutar el modelo de optimización**:

```bash
python main.py
```

   **Nota**: El script ejecutará el solver con presupuesto de 1000 kEUR y límite de tiempo de 600 segundos (10 minutos).

3. **Resultados**:

   La solución se guarda en `5_results/solutions/solution_optimal.csv` y se genera una auditoría de costos en la consola.

## 📊 Flujo de Trabajo

1. **Preparación de Datos** (`data_preparation.py`):
   - Carga datos geográficos desde `2_data/raw/dataset.geojson`
   - Calcula puntuaciones de idoneidad para cada especie según tipo de cobertura del suelo
   - Identifica vecinos de cada celda
   - Exporta dataset final en CSV y GeoJSON

2. **Modelo de Optimización** (`optimization_model.py`):
   - Construye grafo de conectividad entre celdas
   - Calcula rutas mínimas (Dijkstra) desde poblaciones existentes
   - Resuelve problema de optimización con CP-SAT (Google OR-Tools)
   - Considera restricciones de presupuesto, equidad biológica y conflictos entre especies

3. **Análisis y Visualización** (`3_notebooks/`):
   - Análisis exploratorio de datos (EDA)
   - Visualización de mapas y resultados
   - Comparación de iteraciones del modelo

## 🔧 Configuración del Modelo

### Parámetros Principales (en `optimization_model.py`)

- **Presupuesto**: Definido en `main.py` (por defecto: 1000 kEUR)
- **Pesos de especies**: `W_VALS = [1.0, 1.0, 2.0, 1.5]`
- **Equidad biológica**: Rangos mínimos y máximos de área activa por especie
- **Penalización por estrés**: Penaliza conflictos entre especies (Martes + Oryctolagus)

### Restricciones del Modelo

1. **Presupuesto**: Costos totales (adaptación + corredores) ≤ presupuesto
2. **Equidad**: Cada especie debe tener entre X% y Y% del área total activa
3. **Conectividad**: Las celdas activas deben estar conectadas a poblaciones existentes
4. **Conflictos**: Martes y Eliomys no pueden coexistir en la misma celda

## 📈 Resultados

Los resultados incluyen:

- **Solución óptima**: CSV con celdas seleccionadas y especies activas
- **Mapas de visualización**: Expansión de hábitat por especie
- **Auditoría de costos**: Desglose de gastos en adaptación y corredores
- **Métricas de rendimiento**: Score objetivo, gap de optimalidad, tiempo de ejecución

## 📚 Documentación Adicional

Para más detalles técnicos, consulta:

- `DOCUMENTATION.md`: Documentación técnica completa
- `1_documentation/model-design/`: Diseño y evolución del modelo
- `3_notebooks/01_eda_analysis.ipynb`: Análisis exploratorio detallado

## 👥 Autor

Proyecto académico de conservación de biodiversidad en Menorca.
