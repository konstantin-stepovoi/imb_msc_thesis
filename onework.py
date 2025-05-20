import os
import argparse
from JFL3 import one_chromosome_rework_new, process_bed_files  # Импортируем вашу функцию из JFL.py

def process_files(name1, name2, iterations = 50) -> None:
    """
    Обрабатывает одну пару файлов с заданным количеством итераций.
    """
    q = int(name1.split('.')[1])
    try:
        # Получаем Domens и Genes с помощью функции из JFL.py
        Domens, Genes, expressions = process_bed_files(name2, name1)
        print('Splitted files')
        print(f'Domens: {len(Domens)}, {type(Domens), type(Domens[0])}')
        print(f'Genes: {len(Genes)}, {type(Genes), type(Genes[0])}')
        print(f'expresions: {len(expressions)}, {type(expressions), type(expressions[0])}')

        # Запускаем основную обработку
        with open(f"milliontest_1K/{q:05d}.txt", "a") as result_file:
            result_file.write("z_crit\tmu_orig\td_orig\tmu_hor\td_hor\tmu_vert\td_vert\n")
            for _ in range(30): #кол - во повторов задаётся тута
                # Обороты задать не забудь!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                results = one_chromosome_rework_new(Domens, Genes, expressions, 50) #Обороты задаются вот тута!!!
                result_file.write("\t".join(map(str, results)) + "\n")
        print(f'Processed file pair {q}')
    except Exception as e:
        print(f"Error processing {name1}: {e}")

def main():
    """
    Главная функция для чтения аргументов и запуска обработки.
    """
    parser = argparse.ArgumentParser(description="Process one pair of files.")
    parser.add_argument("genes_file", help="Path to the genes file")
    parser.add_argument("domains_file", help="Path to the domains file")
    parser.add_argument("iterations", type=int, help="Number of iterations")
    args = parser.parse_args()

    # Создаём папку, если её нет
    if not os.path.exists("milliontest_1K"):
        os.makedirs("milliontest_1K")

    # Обрабатываем файлы
    process_files(args.genes_file, args.domains_file, args.iterations)

if __name__ == "__main__":
    main()
