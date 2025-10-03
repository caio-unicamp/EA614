import numpy as np
import matplotlib.pyplot as plt


def coef_fourier(k: int):
    '''
    Função pra retornar o coeficiente de ordem k referente à série de fourier do sinal x(t)
    '''
    alfa = k*np.pi/3
    return (1/3)*((np.sin(alfa)/alfa) + (2*np.sin(2*alfa)/alfa) + (2/(alfa**2))*(np.cos(2*alfa) - np.cos(alfa)))

def serie_truncada(N: int, )