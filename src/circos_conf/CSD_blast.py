# Extract circos links from blast hits that fall into CSD loci.
# Cyril Matthey-Doret
# 24.11.2017

import pandas as pd
c_path = "data/circos/"
tig = c_path + "csd_tig.lf.txt"
blast = c_path + "blast.lf.txt"
csd = c_path + "csd_blast.lf.txt"
out = c_path + "out_blast.lf.txt"

tig = pd.read_csv(tig, sep=" ", header=None)
blast = pd.read_csv(blast, sep=" ", header=None)
filtered = pd.DataFrame()
for n in range(tig.shape[0]):
    filtered = filtered.append(blast.loc[(blast[0] == tig.iloc[n,0]) & \
                         (blast[1] >= tig.iloc[n,1]) & \
                         (blast[2] <= tig.iloc[n,2])])
confirmed = pd.DataFrame()
for n in range(tig.shape[0]):
    confirmed = confirmed.append(filtered.loc[(filtered[3] == tig.iloc[n,0]) & \
                         (filtered[4] >= tig.iloc[n,1]) & \
                         (filtered[5] <= tig.iloc[n,2])])
not_csd = blast.drop(confirmed.index)
confirmed.to_csv(csd, sep=" ", index=False, header=False)
not_csd.to_csv(out, sep=" ", index=False, header=False)
