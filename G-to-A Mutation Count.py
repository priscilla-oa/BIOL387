import pandas as pd
g_to_a_path = "/Users/priscilla/Documents/Dry Lab/G_to_A_mutations.xlsx"
g_to_a_df = pd.read_excel(g_to_a_path)
g_to_a_count = g_to_a_df.groupby('#Symbol').size().rename("G_to_A_Mutations").reset_index()
merged_df.to_excel("/Users/priscilla/Documents/Dry Lab/Merged_Gene_Data.xlsx", index=False)


