import numpy as np
import matplotlib.pyplot as plt

# Constantes
T = 6.0 # Período de 6s
w0 = 2*np.pi/T  # omega_0 = 2pi/T
wc = 10 # wc = 1/RC = 1/10^-1 = 10 rad/s

def coef_fourier(k: int):
    '''
    Função pra retornar o coeficiente de ordem k referente à série de fourier do sinal x(t)
    '''
    if k == 0:
        return 2/3.0

    alfa = k*w0
    return (1.0/3.0)*((np.sin(alfa)/alfa) + (2*np.sin(2*alfa)/alfa) + (2/(alfa**2))*(np.cos(2*alfa) - np.cos(alfa)))

def serie_truncada(N: int, tempo):
    '''
    Função que retorna a série truncada de fourier usando os harmônicos de -N a N
    '''
    serie = np.zeros_like(tempo, dtype=complex)

    for k in range(-N, N+1):
        # ~x_n(t) = *sum(a_k*e^(jkw0t))
        a_k = coef_fourier(k) 
        serie += a_k*np.exp(1j * k * w0 * tempo)
    return serie 

def funcao_real(tempo):
    '''
    Função que retorna o valor real do sinal x(t)
    '''
    x_t = np.zeros_like(tempo)

    x_t[(tempo >= -2) & (tempo < -1)] = -2.0 -2.0 *tempo[(tempo >= -2) & (tempo < -1)] 
    x_t[(tempo >= -1) & (tempo < 1)] = 1.0
    x_t[(tempo >= 1) & (tempo < 2)] = -2.0 + 2.0*tempo[(tempo >= 1) & (tempo < 2)]
    return x_t

def resp_freq(omega):
    '''
    Função para saber a resposta em frequência em função de ômega
    '''
    return omega/(omega - 1j*wc)


# Plotagem dos gráficos do item b
qtd_pontos = 2000
tempo = np.linspace(-3, 3, qtd_pontos, False)
dt = tempo[1] - tempo[0]
x_t = funcao_real(tempo)

possiveis_N = [1, 10, 20, 50]
for N in possiveis_N:
    x_n = serie_truncada(N, tempo)
    plt.figure(figsize=(8,3.2))
    plt.plot(tempo, x_t, label='x(t)', color = 'tab:blue')
    plt.plot(tempo, np.real(x_n), label='Série truncada', color = 'tab:orange')
    plt.title(f"Comparação sinal real com a aproximação da série truncada com N = {N}")
    plt.xlabel('t(s)')
    plt.ylabel('Amplitude')
    plt.xlim(-3,3)
    plt.legend()
    plt.grid(True)
    plt.show

# Potência média do erro Pn do item c
for N in possiveis_N:
    xN = np.real(serie_truncada(N, tempo))
    erro = x_t - xN
    p_n = (1.0/T)*np.sum((erro**2)*dt)  # Aproximação para valores discretos não é necessário integrar então aproxima-se numericamente
    print(f"Potência média do erro com N = {N}: {p_n}")

# Módulo dos coeficientes |ak| do item d
indices = np.arange(-50, 51)
coefs = np.array([coef_fourier(k) for k in indices])
omegas = indices*w0
plt.figure(figsize=(8,3.2))
plt.stem(omegas, np.abs(coefs))
plt.title("Módulo dos coeficientes da série de fourier em função de ω")
plt.xlabel("ω (rad/s)")
plt.ylabel("|ak|")
plt.xlim(-12,12)
plt.grid(True)
plt.show

# H(jw) do item e
# Módulo
omega = np.logspace(-2, 3, 1000)
H = resp_freq(omega)
plt.figure(figsize=(8, 3.2))
plt.plot(omega, np.abs(H))
plt.xscale('log')
plt.title("Módulo da Resposta em frequência")
plt.xlabel("ω (rad/s)")
plt.ylabel("|H(ω)|")
plt.grid(True)
plt.show()
# Fase
plt.figure(figsize=(8, 3.2))
plt.plot(omega, np.angle(H, deg=True))
plt.xscale('log')
plt.xlabel("ω (rad/s)")
plt.ylabel("Fase (°) ")
plt.title("Fase da Resposta em frequência")
plt.grid(True)
plt.show()
