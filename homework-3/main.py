import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import csv

POS_PLOT_PATH = "pos_plot.png"
NEG_PLOT_PATH = "neg_plot.png"
POS_DATA_PATH = "pos_data.csv"
NEG_DATA_PATH = "neg_data.csv"

POINTS = 100

T = sp.Symbol('T', real=True)

LAMBDA = 1e-5
beta = 0.0065
betas = [0.0002145, 0.0014235, 0.001274, 0.0025675, 0.0007475, 0.000273]
lambda_vals = [0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01]
lambda1 = 0.0769230


def save_data_to_csv(filename, x_data, y_data_1, y_data_2):
	with open(filename, mode='w', newline='') as file:
		writer = csv.writer(file)
		writer.writerow(["rhos", "six_group", "single_group"])  # Header
		for x, y1, y2 in zip(x_data, y_data_1, y_data_2):
			writer.writerow([x, y1, y2])


def get_ass_t_one_group(rho, is_max=True):
	equation = (LAMBDA / T) + beta / (1 + lambda1 * T) - rho
	solutions = sp.solve(equation, T)
	real_solutions = [sol.evalf() for sol in solutions if sol.is_real]
	print(f'For rho/beta = {rho / beta} -> Tass = {- min(real_solutions)}c')
	if is_max:
		return max(real_solutions)
	else:
		return min(real_solutions)


def get_ass_t_six_groups(rho, is_max=True):
	equation = (LAMBDA / T) + sum(b / (1 + l * T) for b, l in zip(betas, lambda_vals)) - rho
	solutions = sp.solve(equation, T)
	real_solutions = [sol.evalf() for sol in solutions if sol.is_real]
	print(f'For rho/beta = {rho / beta} -> Tass = {- min(real_solutions)}c')
	if is_max:
		return max(real_solutions)
	else:
		return min(real_solutions)


if __name__ == "__main__":
	rhos = np.linspace(0.1, 3, POINTS) * beta
	T_ass_6 = np.array([get_ass_t_six_groups(rho) for rho in rhos])
	T_ass_1 = np.array([get_ass_t_one_group(rho) for rho in rhos])
	save_data_to_csv(POS_DATA_PATH, rhos, T_ass_6, T_ass_1)

	neg_rhos = - np.logspace(-1, 3, POINTS) * beta
	neg_T_ass_6 = np.array([get_ass_t_six_groups(rho, False) for rho in neg_rhos])
	neg_T_ass_1 = np.array([get_ass_t_one_group(rho, False) for rho in neg_rhos])
	save_data_to_csv(NEG_DATA_PATH, neg_rhos, neg_T_ass_6, neg_T_ass_1)

	plt.figure(figsize=(8, 5))
	plt.plot(rhos / beta, T_ass_6, linestyle='-', label="six_group")
	plt.plot(rhos / beta, T_ass_1, linestyle='-', label="single_group")
	plt.yscale("log")  # Логарифмическая шкала только для Y
	plt.xlabel(r"$\rho / \beta$")
	plt.ylabel(r"$T_{ass}$")
	plt.legend()
	plt.grid(True, which="both", linestyle="--")
	plt.savefig(POS_PLOT_PATH)
	plt.show()

	plt.figure(figsize=(8, 5))
	plt.plot(- neg_rhos / beta, - neg_T_ass_6, linestyle='-', label="six_group")
	plt.plot(- neg_rhos / beta, - neg_T_ass_1, linestyle='-', label="single_group")
	plt.yscale("linear")  # Логарифмическая шкала только для Y
	plt.xscale("log")
	plt.xlabel(r"-$\rho / \beta$")
	plt.ylabel(r"-$T_{ass}$")
	plt.legend()
	plt.grid(True, which="both", linestyle="--")
	plt.savefig(NEG_PLOT_PATH)
	plt.show()
