import pandas as pd

def gen_VEP_medium_data_UKBIOBANK(UKBiobank_GWAS_data):

    # filter for SNP
    def if_SNP_bothsex(s):
        return len(s.split(":")[2])==1 and len(s.split(":")[3])==1 and s.split(":")[0]!="X" and s.split(":")[2]!="-" and s.split(":")[3]!="-"
    UKBiobank_GWAS_data = UKBiobank_GWAS_data[UKBiobank_GWAS_data["variant"].apply(if_SNP_bothsex)]

    result = pd.DataFrame() #result VEP metadata

    #generate chromosome, start and allele
    result[["chromosome", "start", "allele"]] = UKBiobank_GWAS_data["variant"].str.split(":",  expand=True, n=2)

    result["start"] = result["start"].astype(int) # generate start
    result["strand"] = "+"                        # generate strand

    #generate end
    # def gen_ends(row):
    #     return row["start"]+len(row["allele"].split(":")[0])-1
    # result["end"] = result.apply(gen_ends, axis=1)
    result["end"] = result["start"]

    #change columns' order
    result = result[["chromosome", "start", "end", "allele", "strand"]]
    return result

def gen_permut_data_UKBIOBANK(UKBiobank_GWAS_data, VEPoutput, p_threshold):
    # generate variants in UKbiobank data, used to match VEP output's first column
    UKBiobank_GWAS_data["variants"] = UKBiobank_GWAS_data["variant"].str.replace(r'([^:]+):([^:]+):([^:]+):([^:]+)', r'\1_\2_\3:\4', regex=True)
    UKBiobank_GWAS_data.drop(columns=["minor_allele", "minor_AF", "expected_case_minor_AC", "low_confidence_variant", "n_complete_samples", "AC", "ytx", "beta", "se", "tstat"], inplace=True)

    # find all annotations in form Consequence_IMPACT， and create new columns for Consequence_IMPACT
    VEPoutput.drop_duplicates(inplace=True)
    VEPoutput.drop(columns=['Location', 'Allele'], inplace=True)
    VEPoutput['Consequence'] = VEPoutput['Consequence'].str.split(',')
    VEPoutput = VEPoutput.explode("Consequence")
    VEPoutput["annotation"] = VEPoutput["Consequence"]+"_"+VEPoutput["IMPACT"]
    VEPoutput.drop(columns=["Consequence", "IMPACT"], inplace=True)
    dummies = pd.get_dummies(VEPoutput['annotation'], prefix='', prefix_sep='')
    VEPoutput.drop(columns=["annotation"], inplace=True)
    agg_dict = {"Uploaded_variation" :"first"}
    for i in dummies.columns:
        agg_dict[i] = "max"
    VEPoutput = pd.concat([VEPoutput, dummies], axis=1).groupby("Uploaded_variation", sort=False).agg(agg_dict)
    VEPoutput.reset_index(drop=True, inplace=True)

    # merge GWAS data and annotation data
    res = UKBiobank_GWAS_data.merge(VEPoutput, left_on='variants', right_on='Uploaded_variation', how="right")
    res.drop(columns=["variants", "Uploaded_variation"], inplace=True)
    return res