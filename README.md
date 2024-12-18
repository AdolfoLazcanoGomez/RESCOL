# Ant Colony Optimization para recolección de residuos domiciliarios RESCOL

## Instrucciones

Compilar en sistemas UNIX con:

`make clean`

`cmake .`

`make`

Ejecucion basica con:

`./RESCOL Instancias/Formato5x5.txt --metodo 0 --eval-max 3000 --num-hormigas 20 --epocas 1 --matriz-secundaria --beta0 --usar-evaluaciones --usar-limitador --valor-limitador 2 --silence --random-seed --dir-salida ./resultados --prefijo-salida nuevos --rescol `


Ejecutar con 10 segundos configurados:

`./RESCOL Instancias/Formato5x5.txt --metodo 0 --num-hormigas 4 --salida-dijkstra --beta0 --alfa 4.29 --rho 0.29 --tau-as 5.07 --valor-limitador 4 --epocas 1 --valor-sin-nuevas-visitas 1 --usar-tiempo --tiempo-max 10  --silence false --dir-salida ./resultados --prefijo-salida nuevos --rescol`

Ejecutar con 30 segundos configurados:

`./RESCOL Instancias/Formato5x5.txt --metodo 0 --num-hormigas 23 --salida-dijkstra --beta0 --alfa 4.11 --rho 0.3 --tau-as 3.57 --valor-limitador 2 --epocas 1 --valor-sin-nuevas-visitas 1 --usar-tiempo --tiempo-max 30  --silence false --dir-salida ./resultados --prefijo-salida nuevos --rescol`

Ejecutar con 1 minuto configurado:

`./RESCOL Instancias/Formato5x5.txt --metodo 0 --num-hormigas 16 --salida-dijkstra --beta0 --alfa 2.79 --rho 0.29 --tau-as 2.39 --valor-limitador 3 --epocas 1 --valor-sin-nuevas-visitas 1 --usar-tiempo --tiempo-max 60  --silence false --dir-salida ./resultados --prefijo-salida nuevos --rescol`


Ejecutar con 5 minutos configurado:
`./RESCOL Instancias/Formato5x5.txt --metodo 0 --num-hormigas 26 --salida-dijkstra --beta0 --alfa 1.81 --rho 0.22 --tau-as 2.14 --valor-limitador 3 --epocas 1 --valor-sin-nuevas-visitas 1 --usar-tiempo --tiempo-max 300  --silence false --dir-salida ./resultados --prefijo-salida nuevos --rescol`

## Parámetros
### Entrada y salida
* --prefijo-salida : Prefijo para archivos de salida [default: resultados]
* --dir-salida : Directorio para archivos de salida [default: ./resultados]
* --rescol : formato entrada y salida para la interfaz rescol [implicit: "true", default: false]
* --irace : usa el formato de llamada e instancias de irace [implicit: "true", default: false]
* --silence : output minimo [implicit: "true", default: false]
* --bdconn : Flag que permite o no subir los resultados a la BD [implicit: "true", default: false]
### Generador de números aleatorios
* --semilla : semilla de generacion aleatoria [default: 1188444438]
* --random-seed : semilla aleatoria [implicit: "true", default: false]
### Budget
* --iter-max : Cantidad de iteraciones máximas [default: 1000]
* --eval-max : Cantidad de evaluaciones máximas [default: 500]
* --tiempo-max: Cantidad maxima de tiempo en segundos [default: 30]
* --usar-iteraciones : usar iteraciones,no se puede usar al mismo tiempo que las demas [implicit: "true", default: false]
* --usar-evaluaciones : usa evaluaciones, no se puede usar al mismo tiempo que las demas [implicit: "true", default: false]
* --usar-tiempo: usa tiempo, no se puede usar al mismo tiempo que las demas. [implicit: "true", default: false]
* --epocas : Épocas [default: 3]
### Base
* --metodo : Metodo de resolucion [default: 1]
* --num-hormigas : Número de hormigas [default: 10]
* --umbral-inf : Umbral inferior para las feromonas [default: 0.001]
* --tau-as : Parámetro tau, asociado a las feromonas iniciales [default: 1]
* --alfa : Parámetro alfa [default: 1]
* --beta : Parámetro beta [default: 2]
* --beta0 : Flag que fuerza el valor de beta a 0 en el recorrido [implicit: "true", default: false]
* --rho : Parámetro rho, asociado a la evaporación de feromonas [default: 0.5]
* --rho-sec : Parámetro rho, asociado a la evaporación de feromonas [default: 0.8]
### Matriz de salida
* --matriz-secundaria : Flag que permite o no el uso de la matriz de salida [implicit: "true", default: false]
* --beta-salida : Parámetro beta matriz de salida [default: 2]
* --rho-salida : Parámetro rho matriz de salida, asociado a la evaporación de feromonas [default: 0.8]
### Algoritmo de dijkstra
* --salida-dijkstra : Flag que permite el uso del algoritmo de dijkstra para salir de la ruta, no se puede al mismo tiempo que la matriz de salida [defaul: false]
### Limitador de pasadas
* --usar-limitador : usa el limitador de pasadas [implicit: "true", default: false]
* --valor-limitador : cantidad limite de veces que se puede pasar por un arco [default: 4]
### Sin nuevas visitas
* --valor-sin-nuevas-visitas: Numero de pasadas por enlaces sin encontrar un enlace nuevo. Antes de recorrer de forma aleatoria hasta encontrar un enlace nuevo (que es seleccionado de forma obligatoria). [default: 1]
### MMAS
* --umbral-sup : Umbral superior para las feromonas [default: 161161]
* --umbral-sin-mejora : Cantidad de iteraciones sin mejora para la actualización de feromonas [default: 10]
* --a-mm : Parámetro a, asociado a la actualización de feromonas [default: 13]
### EAS
* --factor-elitist: Parametro asociado a la actualizacion de feromonas del algoritmo elitista [default: 1.0]
### ACS
* --q0 : Parámetro q0 [default: 0.389]
* --csi : Parámetro csi [default: 0.1]
### Oscilador
* --oscilador : Parametro que permite o no el uso del oscilador de parametros, 0 es desactivado, 1 es stepper, 2 es caotico [default: 0]
* --inc-alfa : Incremento de alfa [default: 0.1]
* --min-alfa : Valor minimo de alfa [default: 1]
* --max-alfa : Valor maximo de alfa [default: 5]
* --inc-beta : Incremento de beta [default: 0.1]
* --min-beta : Valor minimo de beta [default: 1]
* --max-beta : Valor maximo de beta [default: 0.5]
* --inc-rho : Incremento de rho [default: 0.1]
* --min-rho : Valor minimo de rho [default: 0.1]
* --max-rho : Valor maximo de rho [default: 0.5]
* --inc-rho-sec : Incremento de rho secundario [default: 0.1]
* --min-rho-sec : Valor minimo de rho secundario [default: 0.1]
* --max-rho-sec : Valor maximo de rho secundario [default: 0.5]
