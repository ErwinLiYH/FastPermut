import time
import numpy as np

def separate(chromosome_array:np.ndarray, location_array:np.ndarray, annot_array:np.ndarray):
    unique_chromosomes = np.unique(chromosome_array)
    location_groups = [location_array[chromosome_array == chromosome] for chromosome in unique_chromosomes]
    annot_groups = [annot_array[chromosome_array == chromosome] for chromosome in unique_chromosomes]
    return location_groups, annot_groups

def compare_arrays(arr1:np.ndarray, arr2:np.ndarray):
    start_time = time.time()
    if arr1.shape != arr2.shape:
        raise ValueError("Arrays must have the same shape to be compared")

    comparison = arr1 == arr2
    match_count = np.count_nonzero(comparison)
    total_elements = comparison.size

    match_percentage = (match_count / total_elements) * 100
    print(f'Match Percentage: {match_percentage}%', flush=True)
    end_time = time.time()
    print(f"comparison took {end_time - start_time} seconds to run.", flush=True)
    if match_count == total_elements:
        return True
    else:
        return False

def test(to_test:np.ndarray, RNs:np.ndarray, location_vecs:list, annot_mats:list, roll_direction="R", roll_inplace=False, interval=1000):
    start_time = time.time()
    N_times = RNs.shape[1]
    N_chros = RNs.shape[0]
    N_annot = annot_mats[0].shape[1]
    if len(location_vecs)!=N_chros:
        raise Exception("Incorrect number of chromosomes, %d in RNs but %d in location_vecs"%(N_chros, len(location_vecs)))
    if len(annot_mats)!=N_chros:
        raise Exception("Incorrect number of chromosomes, %d in RNs but %d in annot_mats"%(N_chros, len(annot_mats)))
    py_result = np.zeros((N_chros*N_times, N_annot), dtype=np.float32)
    # begin permutation
    for idx, (l, a) in enumerate(zip(location_vecs, annot_mats)): # loop on chromosome
        for i in range(N_times): # loop of permutation
            #print process
            if i!=-1:
                if i%interval==0:
                    print("for chromosome %d, permutation %d"%(idx,i), flush=True)

            # get random value in Stored Random Values
            r = RNs[idx, i]

            if roll_direction=="R":# keep r if roll to right
                pass
            elif roll_direction=="L":# adjust r if roll to left
                r = len(l)-r

            #get rolled vector
            new_l = np.roll(l, r)

            if roll_inplace==False:# keep l(location vector)
                pass
            else:# replace l after geting rolled vector
                l = new_l

            # store result of this permutation to py_result
            py_result[idx*N_times+i, :] = np.matmul(a.T, new_l)

    # compare the result
    compare_arrays(py_result, to_test)
    end_time = time.time()
    print(f"correctness tesing took {end_time - start_time} seconds to run.", flush=True)
    return py_result