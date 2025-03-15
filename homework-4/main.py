import os
import pandas as pd


def read_txt_file(file_path):
	"""Читает .txt файл и возвращает DataFrame с двумя столбцами."""
	data = []
	with open(file_path, 'r', encoding='utf-8') as file:
		for line in file:
			parts = line.split()
			if len(parts) == 2:
				try:
					data.append((float(parts[0]), float(parts[1])))
				except ValueError:
					continue  # Пропускаем некорректные строки

	return pd.DataFrame(data, columns=['T', 'X'])


def merge_txt_files_to_excel(output_file='merged_data.xlsx'):
	"""Ищет все .txt файлы в папке, объединяет их данные в Excel, каждый файл в новые столбцы."""
	all_data = []
	file_names = []

	for file_name in os.listdir('.'):  # Ищем файлы в текущей директории
		if file_name.endswith('.txt'):
			df = read_txt_file(file_name)
			if not df.empty:
				df.columns = [f'T_{file_name}', f'X_{file_name}']  # Добавляем имя файла в название столбцов
				all_data.append(df)
				file_names.append(file_name)

	if all_data:
		merged_df = pd.concat(all_data, axis=1)  # Объединяем данные по ширине, не обрезая строки

		with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
			merged_df.to_excel(writer, index=False)

		print(f'Файл {output_file} успешно создан.')
	else:
		print('Не найдено подходящих .txt файлов.')


if __name__ == "__main__":
	merge_txt_files_to_excel()