import numpy as np
from scipy.special import spherical_jn, spherical_yn
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET


class EPRCalculator:
    def __init__(self, D: float, c=3e8, N: int = 200):
        self.c = c
        self.N = N
        self.r = D / 2

    def bessel(self, n, k, r):
        return spherical_jn(n, k * r) + 1j * spherical_yn(n, k * r)

    def an(self, n, k, r):
        return spherical_jn(n, k * r) / self.bessel(n, k, r)

    def bn(self, n, k, r):
        numerator = k * r * spherical_jn(n - 1, k * r) - n * spherical_jn(n, k * r)
        denominator = k * r * self.bessel(n - 1, k, r) - n * self.bessel(n, k, r)
        return numerator / denominator

    def sigma(self, lmbd, k, r):
        s = 0
        for n in range(1, self.N + 1):
            s += ((-1) ** n * (n + 0.5) * (self.bn(n, k, r) - self.an(n, k, r)))
        return (lmbd ** 2 / np.pi) * abs(s) ** 2

    def calculate(self, frequencies) -> list:
        sigmas = []
        for f in frequencies:
            lmbd = self.c / f
            k = 2 * np.pi / lmbd
            sigmas.append(self.sigma(lmbd, k, self.r))
        return sigmas


def save_to_txt(frequencies, sigmas, filename="result.txt"):
    with open(filename, "w", encoding="utf-8") as f:
        for freq, rcs in zip(frequencies, sigmas):
            freq_str = f"{freq:.6e}".replace('e+0', 'e+').replace('e-0', 'e-')
            rcs_str = f"{rcs:.6e}".replace('e+0', 'e+').replace('e-0', 'e-')
            f.write(f"{freq_str}    {rcs_str}\n")


if __name__ == '__main__':
    input_xml = 'task_rcs_02.xml'
    output_txt = 'result.txt'
    variant_numb = 9

    tree = ET.parse(input_xml)
    root = tree.getroot()

    for variant in root.findall('variant'):
        if int(variant.get('number')) == variant_numb:
            D = float(variant.find('D').text)
            fmin = float(variant.find('fmin').text)
            fmax = float(variant.find('fmax').text)
            break

    print(f"Вариант {variant_numb}: D={D} м, fmin={fmin} Гц, fmax={fmax} Гц")
    frequencies = np.linspace(fmin, fmax, 1000)
    calculator = EPRCalculator(D=D, N=50)
    sigmas = calculator.calculate(frequencies)

    save_to_txt(frequencies, sigmas, filename=output_txt)
    print(f"Результаты сохранены в {output_txt}")

    plt.plot(frequencies / 1e9, sigmas)
    plt.xlabel('Частота (ГГц)')
    plt.ylabel('ЭПР (м²)')
    plt.title(f'ЭПР сферы (D={D} м)')
    plt.grid(True)
    plt.savefig('rcs_epr.png')
    plt.show()
