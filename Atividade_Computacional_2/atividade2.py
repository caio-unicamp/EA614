import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt


def coef_fourier(k: int):
    '''
    Função pra retornar o coeficiente de ordem k referente à série de fourier do sinal x(t)
    '''
    if k == 0:
        return 2/3

    alfa = k*np.pi/3
    return (1/3)*((np.sin(alfa)/alfa) + (2*np.sin(2*alfa)/alfa) + (2/(alfa**2))*(np.cos(2*alfa) - np.cos(alfa)))

def serie_truncada(N: int, tempo):
    '''
    Função que retorna a série truncada de fourier usando os harmônicos de -N a N
    '''
    x_n = 0
    for k in range(-N, N+1):
        x_n += coef_fourier(k)*np.exp(1j*k*np.pi*tempo/3) # ~x_n(t) = e^(t)*sum(a_k*e^(jkw0))

    return x_n

def funcao_real(tempo):
    '''
    Função que retorna o valor real do sinal x(t)
    '''
    if -2<= tempo< 1:
        return -2 -2*tempo
    elif -1 <= tempo < 1:
        return 1
    elif 1<= tempo < 2:
        return -2 + 2*tempo
    else:
        return 0

def erro_Pn(N: int):
    return (1/6)* integrate((funcao_real - serie_truncada(N = N))**2, -3, 3)    # Pn = 1/T int((x(t) - ~x_n(t))²dt) de -T/2 a T/2