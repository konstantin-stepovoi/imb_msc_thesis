from DomenSorter import Domen_Sort_f

a = Domen_Sort_f('Gentype_uni.txt', 0)
b = Domen_Sort_f('Gentype_uni.txt', 1)

with open('domens_uni.txt', 'a') as file:
    file.write(a)
    file.write(b)

print('Done')