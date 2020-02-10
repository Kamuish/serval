from src import builder

import os
import glob

files = glob.glob('gj699/*')


for f in files:
    os.remove(f)

builder()


print("EVERYTHING IS DONE !!!! WE REACHED THE END")