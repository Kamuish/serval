from src import builder
from src import main

import os, sys
import glob


plot = 1

if plot:
    main()
    sys.exit(0)


files = glob.glob('gj699/*')

for f in files:
    os.remove(f)

builder()


print("EVERYTHING IS DONE !!!! WE REACHED THE END")