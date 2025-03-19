import pandas as pd
excel_path = r"C:\Users\prige\OneDrive\Desktop\Dry Lab\Dry Lab Data.xlsx"
data = pd.read_excel(excel_path)
print(data.head())
g_to_a_mutations = data[data['NucleotideChange'].str.contains('G>A', na=False)]
print(g_to_a_mutations)
g_to_a_mutations.to_excel(r"C:\Users\prige\OneDrive\Desktop\Dry Lab\G_to_A_mutations.xlsx", index=False)
import os
print("Current Working Directory:", os.getcwd())
g_to_a_mutations = data[data['NucleotideChange'].str.contains('G>A', na=False)]
file_path = r"C:\Users\prige\OneDrive\Desktop\G_to_A_mutations.xlsx"
g_to_a_mutations.to_excel(file_path, index=False)
print(f"File saved to {file_path}")

#types of mutations
mutation_counts = g_to_a_mutations['mutation_consequence'].value_counts()
print(mutation_counts)




