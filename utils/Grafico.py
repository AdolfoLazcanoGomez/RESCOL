import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sys import argv

# Lee el archivo, asumiendo que el formato es exactamente 'Epoca: %d, Evaluacion: %d, Mejor costo: %d'
#arg = 'Output/RencaChiquito-AntSystem-20231020232543/20231020232543.txt'
df = pd.read_csv(argv[1], sep=", ", engine='python', header=None, names=['Epoca', 'Evaluacion', 'Mejor costo'])

# Convierte las columnas al tipo de datos correcto y extrae los números de las cadenas
df['Epoca'] = df['Epoca'].str.extract('(\d+)').astype(int)
df['Evaluacion'] = df['Evaluacion'].str.extract('(\d+)').astype(int)
df['Mejor costo'] = df['Mejor costo'].str.extract('(\d+)').astype(int)

# Filtra la fila con el valor máximo de 'Mejor costo'
df = df[df['Mejor costo'] != df['Mejor costo'].max()]

# Agrupa por 'Epoca' e 'Evaluacion', tomando el mínimo 'Mejor costo' para cada grupo
df_grouped = df.groupby(['Epoca', 'Evaluacion']).min().reset_index()

# Crea un color distinto para cada época
palette = sns.color_palette('hsv', df_grouped['Epoca'].nunique())

# Crea un gráfico de lineas para cada época
for epoca in df_grouped['Epoca'].unique():
    data = df_grouped[df_grouped['Epoca'] == epoca]
    
    plt.plot(data['Evaluacion'], data['Mejor costo'], label=f'Epoca {epoca}', color=palette[epoca])
    

plt.xlabel('Evaluaciones')
plt.ylabel('Mejor costo')
plt.title('Gráfico de convergencia')
plt.legend()
nombre_archivo = argv[1][:-4]
plt.savefig(nombre_archivo + '-GraficoConvergencia.png')

show_grafico = argv[2]
if (show_grafico == "True"):
    plt.show()
