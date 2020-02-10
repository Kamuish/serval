
import tarfile
from astropy.io import fits 


class FitsClass():
    def __init__(self, filepath, name, offset_data=0, size=0, mode = 1, **kwargs):
        """
            mode:
                1 -> going through the file manually 
                2 -> tarfile stored inside

        """
        try:
            self.file_path = [filepath, kwargs['tar_file']]
        except:
            self.file_path = filepath
        self.open_file = None

        self.name = name 

        self.offset_data = offset_data 

        self.size = size 

        self.mode = mode

    def __enter__(self):
        if self.open_file is not None:
            return self.open_file

        if self.mode == 1:
            self.open_file = open(self.file_path, mode = 'rb')
        elif self.mode == 2:
            self.open_file = fits.open(self.file_path[0])
        return self.open_file 

    def __exit__(self, type, value, traceback):
        if self.mode != 2:   
            self.open_file.close()
            self.open_file = None       


    def __del__(self):
        """
            Makes sure that the file is closed
        """
        if self.open_file is not None and self.mode != 2:
            self.file_path.close()
        if self.mode == 2:
            print("Closing down the tar file open in memory")
            self.file_path[0].close()
            self.file_path[1].close()

    def seek(self, seek_pos, whence = 0):
        self.open_file.seek(int(seek_pos), whence)

   
    
if __name__ == '__main__':
    import tarfile 
    s = '/home/amiguel/SERVALHOME/serval/data/HARPS/gj699/ADP.2014-09-17T11-23-26.740.tar'
    tar = tarfile.open(s)

    for member in tar.getmembers():
      if '_e2ds_A.fits' in member.name: 
        extr = member

        s = FitsClass(filepath = s, 
                offset_data = extr.offset_data,
                name = extr.name,
                size=extr.size) 


    s.seek(s.offset_data)    
    args = (b'INSTRUME', b'OBJECT', b'MJD-OBS')
    for card in iter(lambda:s.file.read(80), ''):   # read in 80 byte blocks
        if card.startswith(args):
            print('yes')

    tar.close()