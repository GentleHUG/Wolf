import os
import pandas as pd


def read_txt_file(file_path):
	"""Читает .txt файл и возвращает DataFrame с двумя столбцами."""
	data = []
	with open(file_path, 'r', encoding='utf-8') as file:
		for line in file:
			parts = line.split()
			if len(parts) == 2:  # Убеждаемся, что в строке два значения
				try:
					data.append((float(parts[0]), float(parts[1])))
				except ValueError:
					continue  # Пропускаем строки, которые не удаётся конвертировать

	return pd.DataFrame(data, columns=['X', 'Y'])


def merge_txt_files_to_csv(output_file='merged_data.csv'):
	"""Ищет все .txt файлы в папке, объединяет их данные и сохраняет в CSV."""
	all_data = []

	for file_name in os.listdir('.'):  # Ищем файлы в текущей директории
		if file_name.endswith('.txt'):
			df = read_txt_file(file_name)
			if not df.empty:
				df['File'] = file_name  # Добавляем имя файла как дополнительный столбец
				all_data.append(df)

	if all_data:
		merged_df = pd.concat(all_data, ignore_index=True)
		merged_df.to_csv(output_file, index=False, encoding='utf-8')
		print(f'Файл {output_file} успешно создан.')
	else:
		print('Не найдено подходящих .txt файлов.')


if __name__ == "__main__":
	merge_txt_files_to_csv()
