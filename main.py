import pandas as pd
import re

BETA_FILE_PATH = "beta.txt"
NEUT_FILE_PATH = "neut.txt"

# Разбиение на отрезки, по которым будет произодиться агрегация
BINS = [0, 400, 2_000, 5_000, 20_000, 80_000]
labels = [f'{BINS[i]}-{BINS[i + 1]}' for i in range(len(BINS) - 1)]


def parse_beta_file(file_path):
	data = []
	with open(file_path, 'r', encoding='utf-8') as file:
		for line in file:
			parts = re.split(r'\s+', line.strip())
			if len(parts) == 6:
				isotope, value1, value2, unit, value3, value4 = parts
				data.append([isotope, float(value3), float(value4)])

	df = pd.DataFrame(data, columns=["Isotope", "HalflifeMs", "NeutPerDecay"])
	return df[df["NeutPerDecay"] != 0]


def parse_neut_file(file_path):
	data = []
	with open(file_path, 'r', encoding='utf-8') as file:
		for line in file:
			parts = re.split(r'\s+', line.strip())
			if len(parts) == 3:
				isotope, value1, value2 = parts
				data.append([isotope, float(value1)])

	df = pd.DataFrame(data, columns=["Isotope", "Frequency"])
	return df


def group_isomeric_isotopes(df):
	df["BaseIsotope"] = df["Isotope"].str.replace(r'(\d+)m$', r'\1', regex=True)
	grouped_df = df.groupby("BaseIsotope").agg({
		"Frequency": "sum"
	}).reset_index()
	grouped_df.rename(columns={"BaseIsotope": "Isotope"}, inplace=True)
	return grouped_df


if __name__ == "__main__":
	beta_df = parse_beta_file(BETA_FILE_PATH)
	print(beta_df)
	neut_df = parse_neut_file(NEUT_FILE_PATH)
	print(neut_df)
	neut_df = group_isomeric_isotopes(neut_df)
	print(neut_df)

	merged_df = pd.merge(beta_df, neut_df, on="Isotope", how="inner")
	merged_df = merged_df.sort_values(by="HalflifeMs", ascending=False)

	merged_df['HalflifeMsGroup'] = pd.cut(merged_df['HalflifeMs'], bins=BINS, labels=labels, right=False)

	print(merged_df)

	merged_df['Frequency'] = merged_df['Frequency'] / merged_df['Frequency'].sum()
	merged_df.rename(columns={"Frequency": "WeightedFrequency"}, inplace=True)

	merged_df['beta'] = merged_df["WeightedFrequency"] * merged_df["NeutPerDecay"]


	grouped = merged_df.groupby('HalflifeMsGroup', observed=False)

	def weighted_avg(group):
		weighted_halflife = (group['HalflifeMs'] * group['WeightedFrequency']).sum() / group['WeightedFrequency'].sum()
		group_beta = group['beta'].sum()
		return pd.Series({"WeightedAvgHalflifeMs": weighted_halflife, "GroupBeta": group_beta})

	grouped_data = grouped[['HalflifeMs', 'WeightedFrequency', 'beta']]

	groups = grouped_data.apply(weighted_avg)

	groups['WeightedBeta'] = groups['GroupBeta'] / groups['GroupBeta'].sum()

	print(groups[['WeightedAvgHalflifeMs', 'WeightedBeta']])
	print(groups['GroupBeta'].sum())