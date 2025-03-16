import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import xlabel


def load_data(filename):
	data = np.loadtxt(filename, skiprows=1)  # Пропускаем первую строку
	return data[:, 0], data[:, 1]


def plot_data(x, y, output_filename, title, xlabel, ylabel):
	plt.figure(figsize=(10, 5))
	plt.plot(x, y)

	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.grid(True)
	plt.savefig(output_filename)
	plt.show()


def main():
	file1 = 'neut1.txt'  # Замените на ваш файл
	file2 = 'neut2.txt'  # Замените на ваш файл
	neut_plot = 'neut_plot.png'

	title_neut = 'Зависимость плоности нейтронов от времени'
	xlabel_neut = 'T, s'
	ylabel_neut = 'neut'

	t1, neut1 = load_data(file1)
	t2, neut2 = load_data(file2)
	t2 = t2 + t1[-1]

	# Объединяем данные
	x_combined = np.concatenate((t1, t2))
	y_combined = np.concatenate((neut1, neut2))

	plot_data(x_combined, y_combined, neut_plot, title_neut, xlabel_neut, ylabel_neut)

	file1 = 'rho1.txt'  # Замените на ваш файл
	file2 = 'rho2.txt'  # Замените на ваш файл
	neut_plot = 'rho_plot.png'

	title_neut = 'Зависимость реактивности от времени'
	xlabel_neut = 'T, s'
	ylabel_neut = r'$\rho / \beta$'

	t1, neut1 = load_data(file1)
	t2, neut2 = load_data(file2)
	t2 = t2 + t1[-1]

	# Объединяем данные
	x_combined = np.concatenate((t1, t2))
	y_combined = np.concatenate((neut1, neut2))

	plot_data(x_combined, y_combined, neut_plot, title_neut, xlabel_neut, ylabel_neut)



if __name__ == "__main__":
	main()
