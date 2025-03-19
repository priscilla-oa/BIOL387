import pandas as pd
import matplotlib.pyplot as plt


g_to_a_path = "/Users/priscilla/Documents/Dry Lab/G_to_A_mutations.xlsx"
g_to_a_df = pd.read_excel(g_to_a_path)
top_15_genes = ['MYH7', 'MYBPC3', 'FLNC', 'DSP', 'TTN', 'TCAP', 'LZTR1', 'TNNI3', 'COL3A1', 'APOB', 'DOLK', 'SCO2', 'FBN1', 'EMD', 'LDLR']
filtered_df = g_to_a_df[g_to_a_df['#Symbol'].isin(top_15_genes)]
mutation_counts = filtered_df.groupby(['#Symbol', 'mutation_consequence']).size().unstack(fill_value=0)
mutation_counts = mutation_counts.apply(lambda row: row[row != row.max()], axis=1)
mutation_proportions = mutation_counts.div(mutation_counts.sum(axis=1), axis=0)
top_7_genes = top_15_genes[:7]
top_8_genes = top_15_genes[7:]
top_7_proportions = mutation_proportions.loc[top_7_genes]
top_8_proportions = mutation_proportions.loc[top_8_genes]

#Creation of Stacked Bar Chart
fig, ax = plt.subplots(figsize=(10, 6))
top_7_proportions.plot(kind='bar', stacked=True, ax=ax, colormap='tab20', width=0.8)
ax.set_title("Proportions of Mutation Types for Top 7 Genes")
ax.set_xlabel("Gene")
ax.set_ylabel("Proportion of Mutation Types")
ax.legend(title="Mutation Type", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=45)
plt.tight_layout()
fig, ax = plt.subplots(figsize=(10, 6))
top_8_proportions.plot(kind='bar', stacked=True, ax=ax, colormap='tab20', width=0.8)
ax.set_title("Proportions of Mutation Types for Top 8 Genes")
ax.set_xlabel("Gene")
ax.set_ylabel("Proportion of Mutation Types")
ax.legend(title="Mutation Type", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

