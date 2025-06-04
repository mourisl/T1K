import sys
import pandas as pd

# -------------------------------
# Step 1: Read and Insert Intron Rows
# -------------------------------

# Read the initial CSV file containing CFTR transcript data
if len(sys.argv) < 2:
    print("Usage: python3 cftr_exon_intron_coordinate_mapper.py <input_file>")
    sys.exit(1)

input_file = sys.argv[1]
df = pd.read_csv(input_file)

# Determine the last exon number
last_exon_num = df.iloc[-1]['exon_num']

# Create a list to store both the original rows and the new intron rows.
new_rows = []

# Loop through the rows; after each exon row (except the last), insert an intron row.
for i in range(len(df)):
    # Append the current row (as a dict)
    new_rows.append(df.iloc[i].to_dict())
    
    # If the current row is an exon and not the last exon, insert an intron row.
    if df.iloc[i]['type'] == 'exon' and df.iloc[i]['exon_num'] != last_exon_num:
        intron_pos1 = df.iloc[i]['pos2'] + 1
        if i + 1 < len(df):
            intron_pos2 = df.iloc[i + 1]['pos1'] - 1
        else:
            intron_pos2 = None
        intron_row = {
            'chr7': 'chr7',
            'type': 'intron',
            'pos1': intron_pos1,
            'pos2': intron_pos2,
            'exon_num': df.iloc[i]['exon_num'],
            'trans_name': 'CFTR-201'
        }
        new_rows.append(intron_row)

# Create a new DataFrame with both exons and inserted intron rows.
df2 = pd.DataFrame(new_rows)

# Add the bp_num column (number of base pairs) calculated as pos2 - pos1 + 1.
df2['bp_num'] = df2['pos2'] - df2['pos1'] + 1

# -------------------------------
# Add pos1_CFTR and pos2_CFTR Columns
# -------------------------------
# Starting from the first row, set:
#   pos1_CFTR = 1,
#   pos2_CFTR = pos1_CFTR + bp_num - 1.
# Then, for every subsequent row:
#   pos1_CFTR[i] = pos2_CFTR[i-1] + 1,
#   pos2_CFTR[i] = pos1_CFTR[i] + bp_num - 1.

# Initialize lists to store CFTR positions.
pos1_CFTR = []
pos2_CFTR = []

# Iterate over rows in the DataFrame in order.
for idx, row in df2.iterrows():
    if idx == 0:
        current_pos1 = 1
    else:
        current_pos1 = pos2_CFTR[-1] + 1
    current_pos2 = current_pos1 + row['bp_num'] - 1
    pos1_CFTR.append(current_pos1)
    pos2_CFTR.append(current_pos2)

# Add the new columns to the DataFrame.
df2['pos1_CFTR'] = pos1_CFTR
df2['pos2_CFTR'] = pos2_CFTR



# -------------------------------
# Step 2: Define Exon/ Intron Identifiers
# -------------------------------

df2['exon_num'] = df2['exon_num'].fillna(0).astype(int)
# Create a new column "concat" that combines the region type with its exon number (e.g., "exon5" or "intron5")
df2['concat'] = df2['type'] + df2['exon_num'].astype(str)

# Create lists of exons and introns based on the "concat" column.å=
exons = []
introns = []
for _, row in df2.iterrows():
    if row['type'] == 'exon':
        exons.append(row['concat'])
    elif row['type'] == 'intron':
        introns.append(row['concat'])

print("Exons:", exons)
print("Introns:", introns)



# -------------------------------
# Step 3: Map DNA Positions to RNA and cDNA Coordinates
# -------------------------------

# Ensure these DNA positions are integers.
df2['pos1_CFTR'] = df2['pos1_CFTR'].fillna(0).astype(int)
df2['pos2_CFTR'] = df2['pos2_CFTR'].fillna(0).astype(int)

# Initialize RNA position columns.
df2['pos1_RNA'] = None
df2['pos2_RNA'] = None

# Map DNA positions to RNA positions for exons only.
# For the first exon, RNA positions equal the DNA positions.
# For subsequent exons, RNA positions are contiguous, starting one bp after the previous exon's RNA end.
for index, row in df2.iterrows():
    if row['concat'] in exons:  # Process only exon rows
        exon_index = exons.index(row['concat'])
        if exon_index == 0:
            df2.loc[index, 'pos1_RNA'] = row['pos1_CFTR']
            df2.loc[index, 'pos2_RNA'] = row['pos2_CFTR']
        else:
            # Get the previous exon row (using the "concat" identifier)
            previous_exon = exons[exon_index - 1]
            previous_row = df2[df2['concat'] == previous_exon].iloc[0]
            pos1_RNA = previous_row['pos2_RNA'] + 1
            pos2_RNA = pos1_RNA + row['bp_num'] - 1  # Use bp_num from current exon
            df2.loc[index, 'pos1_RNA'] = pos1_RNA
            df2.loc[index, 'pos2_RNA'] = pos2_RNA

# Convert RNA positions to a nullable integer type.
df2['pos1_RNA'] = df2['pos1_RNA'].astype('Int64')
df2['pos2_RNA'] = df2['pos2_RNA'].astype('Int64')

# Map RNA positions to cDNA variant positions (subtracting 70 bp)
df['pos1_cDNA_name'] = None
df['pos2_cDNA_name'] = None
df2['pos1_cDNA_name'] = df2['pos1_RNA'] - 70
df2['pos2_cDNA_name'] = df2['pos2_RNA'] - 70

# Map DNA positions to Python (0-indexed) positions.
df['pos1_CFTR_py'] = None
df['pos2_CFTR_py'] = None
df2['pos1_CFTR_py'] = df2['pos1_CFTR'] - 1
df2['pos2_CFTR_py'] = df2['pos2_CFTR'] - 1



# -------------------------------
# Step 4: Save Final Output
# -------------------------------

cols_to_keep = ['type', 'exon_num', 'concat', 'bp_num', 'pos1_RNA', 'pos2_RNA', 'pos1_CFTR_py', 'pos2_CFTR_py']
df2 = df2[cols_to_keep].copy()
df2.rename(columns={'exon_num': 'num', 'concat': 'exon-intron_structure'}, inplace=True)
output_file = "CFTR-201_Exon_Intron_Complete_Coordinates.csv"
df2.to_csv(output_file, index=False)
