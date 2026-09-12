import pandas as pd
import numpy as np

datos = pd.read_csv('Bases/anime.csv')
generos = pd.read_csv('Bases/anime.csv', names=['anime_id','members','type','episodes','rating'])

# Creamos la matriz de ratings donde:
# - index (filas) = 'members' (usuarios/miembros)
# - columns = 'anime_id' (id del anime)
# - values = 'rating' (calificación del usuario i en anime j)
# - fill_value=0  llena los espacios vacíos con 0
matriz_ratings = datos.pivot_table(
    index='members', 
    columns='anime_id', 
    values='rating', 
    #fill_value=0   
)

# Proporción de datos no observados
total_nans = matriz_ratings.isna().sum().sum()
total_celdas = matriz_ratings.size
proporcion_nans = total_nans / total_celdas
print(f"Proporción de valores no observados: {proporcion_nans:.4f}")
print(f"Es decir, aproximadamente un {proporcion_nans * 100:.2f}% de la matriz está vacío.")

# Guardar la matriz en formato parquet
matriz_ratings.to_parquet('Bases/anime_ratings.parquet')

