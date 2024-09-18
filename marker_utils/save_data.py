from pathlib import Path
import pandas as pd
def save_data(info,output_dir="./output"):
    outdir=Path(output_dir)
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    clusterlist = info.keys()
    for group in clusterlist:
        df=pd.DataFrame(info[group])
        df = df.sort_values(by='EI', ascending=False)
        df.to_csv(f"{outdir}/{group}_EI.csv", index=True)