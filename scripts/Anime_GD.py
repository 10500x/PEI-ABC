import pandas as pd
import numpy as np

datos = pd.read_csv('Bases/anime.csv')
generos = pd.read_csv('Bases/anime.csv', names=['anime_id','members','type','episodes','rating'])


matrizdatos = pd.crosstab(datos['name'], datos['members']).clip(upper=1)
matrizgeneros = pd.crosstab(datos['name'], datos['genre']).clip(upper=1)

matrizdatos.to_parquet('anime_usuario.parquet')
matrizgeneros.to_parquet('anime_genero.parquet')