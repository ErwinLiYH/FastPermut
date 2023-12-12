import numpy as np
import pandas as pd
import os
import glob


# def load_csv_data(data_directory: str, annot_mat_row_or_column_based: str, encoding="utf-8"):
#     chromosome = np.genfromtxt(os.path.join(data_directory, "chromosome.csv"), delimiter=",", skip_header=1,
#                                dtype=np.single, encoding="utf-8")
#     location = np.genfromtxt(os.path.join(data_directory, "location.csv"), delimiter=",", skip_header=1,
#                              dtype=np.single, encoding="utf-8")
#     if annot_mat_row_or_column_based == "R":
#         # numpy.genfromtxt generate row-based array
#         annotation = np.genfromtxt(os.path.join(data_directory, "annotation.csv"), delimiter=",", skip_header=1,
#                                    dtype=np.single, encoding="utf-8")
#     elif annot_mat_row_or_column_based == "C":
#         # pandas.read_csv generate conlumn-based array
#         annotation = pd.read_csv(os.path.join(data_directory, "annotation.csv"), dtype=np.single,
#                                  encoding="utf-8").to_numpy()
#     else:
#         raise Exception('annot_mat_row_or_column_based should be "R"(row) or "C"(column)')
#     return chromosome, location, annotation


def load_from_file(data_path: str, encoding="utf-8", order="F", delimiter=",", glob_str="*.csv"):
    if os.path.isdir(data_path):
        csv_files = glob.glob(os.path.join(data_path, glob_str))
        dfs = []
        for i, j in enumerate(csv_files):
            df = pd.read_csv(j, dtype=np.single, encoding=encoding, delimiter=delimiter)
            df["chromosome"] = i+1
            dfs.append(df)
        combined_df = pd.concat(dfs)
        return combined_df.iloc[:,-1].to_numpy(), combined_df.iloc[:,1].to_numpy(), np.asarray(combined_df.iloc[:,2:-1].to_numpy(), order=order)
    elif os.path.isfile(data_path):
        df = pd.read_csv(j, dtype=np.single, encoding=encoding, delimiter=delimiter)
        return df.iloc[:,0].to_numpy(), df.iloc[:,1].to_numpy(), np.asarray(df.iloc[:,2:-1].to_numpy(), order=order)

def load_RNs(path:str, N_chros:int, N_times:int):
    data = np.fromfile(path, dtype=np.int32)
    if data.shape[0]!=(N_chros*N_times):
        raise Exception("size incorrect! size:%d"%data.shape[0])
    return data.reshape(-1, N_times)

def load_mpi_result(path:str, name:str, N_chros:int, N_annot:int):
    csv_files = glob.glob(os.path.join(path, '%s_R*.csv'%name))
    dfs = []
    total_times=0
    each_times = []
    for i in csv_files:
        df = pd.read_csv(i)
        N_times = int(df.shape[0]/N_chros)
        if df.shape[0]%N_chros!=0:
            raise Exception("Size error! %d records in %s, but number of chromosomes is %d"%(df.shape[0], i, N_chros))
        print("read data in %s, number of iteration: %d"%(i, N_times))
        total_times+=N_times
        each_times.append(N_times)
        dfs.append(df)
    print("Total permutation times: %d"%total_times)
    res_df = pd.DataFrame(columns=[str(i) for i in range(N_annot)])
    for i in range(N_chros):
        for df, times in zip(dfs, each_times):
            res_df = pd.concat([res_df, df.iloc[i*times:(i+1)*times, :]])
    return res_df

def load_VEP_output(path):
    with open(path, 'r') as f:
        for line in enumerate(f):
            if line.startswith('#') and not line.startswith('##'):
                headers = line.strip("#").strip("\n").split("\t")
                break
    res = pd.read_csv(path, header=None, comment='#', sep="\t")
    res.columns = headers
    return res

def load_VEP_outputs(path_func, num_files):
    dfs = []  # List to hold individual DataFrames

    for i in range(num_files):
        file_name = path_func(i)
        if os.path.exists(file_name):
            dfs.append(load_VEP_output(file_name))
            print(f"Loaded file {file_name}")
        else:
            print(f"Warning: {file_name} does not exist.")

    return pd.concat(dfs, ignore_index=True)

def save_df_to_multiple_tsvs(df, path_func, num_files):
    """
    Splits a DataFrame into a specified number of CSVs.

    Parameters:
    - df: DataFrame to be split.
    - path_func: name format for TSV files. if it is "lambda i: f'path/name_{i}.tsv'", files will be path/name_0.tsv, path/name_0.tsv, ...
    - num_files: Number of TSV files to split the DataFrame into.
    """
    rows_per_csv = len(df) // num_files

    for i in range(num_files):
        start_idx = i * rows_per_csv

        # If it's the last file, include all remaining rows. Otherwise, calculate end index.
        if i == num_files - 1:
            end_idx = None
        else:
            end_idx = (i+1) * rows_per_csv

        chunk = df.iloc[start_idx:end_idx]
        chunk.to_csv(path_func(i), index=False, sep="\t", header=False)
        print(f"save to to file {path_func(i)}")