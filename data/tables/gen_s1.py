import pandas as pd

TC_DISPLAY_NAME = {
     "flower"    : "Flower",
     "leaf"      : "Leaf",
     "root"      : "Root",
     "rosette"   : "Rosette",
     "seed"      : "Seed",
     "seedling1wk": "Seedling (1 Week)",
     "seedling2wk": "Seedling (2 Weeks)",
     "shoot"      : "Shoot",
     "wholeplant" : "Whole plant",
     "chemical"   : "Chemical",
     "development" : "Development",
     "hormone-aba-iaa-ga-br" : "Hormone (ABA, IAA, GA, BR)",
     "hormone-ja-sa-ethylene" : "Hormone (JA, SA, Ethylene)",
     "light"        : "Light",
     "nutrients"    : "Nutrients",
     "stress-light" : "Stress (Light)",
     "stress-other" : "Stress (Other)",
     "stress-pathogen"     : "Stress (Pathogen)",
     "stress-salt-drought" : "Stress (Salt, Drought)",
     "stress-temperature"  : "Stress (Temperature)"
}
def get_tcname(tx):
    return TC_DISPLAY_NAME[tx]

scols = ["FileId", "SeriesId", "SeriesTitle", "SeriesLink", "SampleAttributes", "SampleFile"]
out_cols = ["SeriesId", "SeriesTitle", "SeriesLink", "SampleAttributes", "SampleFile"]
msfx = "Dataset-MasterSheet-Accepted.csv"
msupdfx = "GEO-update-list.csv"
out_file = "S1-Dataset-List.xlsx"

msdfmx = pd.read_csv(msfx, low_memory=False)
print("Master", msdfmx.shape)
msudfmx = pd.read_csv(msupdfx, low_memory=False)
msudfmx['FileId'] = msudfmx["SeriesId"].map(lambda x : str(x) + "_") + msudfmx["SampleId"]
print("Update", msudfmx.shape, sum(1 if x in msudfmx.columns else 0 for x in scols) == 6)

with open("tco.txt") as tcof:
    tcos = [x.strip() for x in tcof.readlines()]
print(tcos)

with pd.ExcelWriter(out_file) as writer:
    for tx in tcos:
        tfx = tx + "/final-list.csv"
        tdx = pd.read_csv(tfx, names=["SeriesId"+tx, "FileId"])
        tdmgx = tdx.merge(msdfmx, on="FileId").loc[:, scols]
        tdux = pd.read_csv(tfx, names=["SeriesId"+tx, "FileId"])
        tdmgux = tdux.merge(msudfmx, on="FileId").loc[:, scols]
        print(tfx, tdx.shape, tdmgx.shape, tdmgux.shape, tdx.shape[0] == tdmgux.shape[0] + tdmgx.shape[0])
        vdf = tdmgx.append(tdmgux)
        vdf = vdf.loc[:, out_cols]
        sname = get_tcname(tx)
        print(sname)
        vdf.to_excel(writer, sheet_name=sname, index=False)

#dx.loc[fx.FileId, ["FileId","SeriesId"]]
#dx.loc[fx.FileId, ["SeriesId", "SeriesTitle", "SeriesLink"]]
#dx.loc[fx.FileId, ["SeriesId", "SeriesTitle", "SeriesLink", "SampleAttributes"]]
#dx.loc[fx.FileId, ["SeriesId", "SeriesTitle", "SeriesLink", "SampleAttributes", "SampleFile"]]
#msx.columns
#fx = pd.read_csv("flower/final-list.csv", names=["Series", "FileId"])
#fx.join(msx).shape
#fx.join(msx, on="FileId").shape
#fx.merge(msx, on="FileId").shape
#
