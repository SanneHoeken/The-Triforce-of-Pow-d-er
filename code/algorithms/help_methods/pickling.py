import pickle
import os

def store_pickle(directory_path, pickled):    
    pickle_path = directory_path +  "/pickle" + str(len(os.listdir(directory_path)))
    pickle_file = open(pickle_path, 'ab')
    pickle.dump(pickled, pickle_file)
    pickle_file.close() 
    return pickle_path


def show_pickle(file_path):
    pickle_file = open(file_path, 'rb')      
    pickled = pickle.load(pickle_file) 
    for keys in pickled: 
        print(keys, '=>', pickled[keys]) 
    pickle_file.close() 


def get_pickle(file_path):
    pickle_file = open(file_path, 'rb')      
    pickled = pickle.load(pickle_file) 
    pickle_file.close() 
    return pickled


def show_all_pickles(directory_path):
    for file in os.listdir(directory_path):
        print(f"#### Pickle from { file } ####")
        show_pickle(directory_path +  "/" + file)
