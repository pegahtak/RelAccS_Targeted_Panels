# this script is to count the median length of different fragments mapped to the same region 
import pandas as pd
import statistics as st


bed_df = pd.read_csv("coad_open.bed", sep="\t")
bed_df = bed_df[["chr", "start", "end"]]
print("\n ---- bed_df head ----")
print(bed_df.head())
frag_df = pd.read_csv("coad_open_frag.csv")
print("\n ---- frag_df head ----")
print(frag_df.head())

print("\n\n Finding medians . . . \n")

median_lst = []
for index1,row1 in bed_df.iterrows():
    size_lst = []
    print("\r {}/{}".format(index1 + 1, len(bed_df)), end="")
    for index2,row2 in frag_df.iterrows():
        if((row2["chr"]==row1["chr"])&(row2["start"]>=row1["start"])&(row2["end"]<=row1["end"])):
            size_lst.append(row2["size"])

    if len(size_lst) > 0:
        median_lst.append(st.median(size_lst))
    else:
        median_lst.append(-1)


bed_df["size_median"] = median_lst
bed_df.to_csv("new_bed_data.csv", index=False)
print("\n\n ---- bed_df head ----")
print(bed_df.head())
